"""
ensemble.py

Tools for generating and writing 3D conformational ensembles with one bead per
residue.

There are two kinds of generator here:

* ``gaussian_chain_factor()`` / ``sample_gaussian_chain()`` work for any
  *Gaussian* chain - one in which every inter-residue distance is
  Gaussian-distributed. Such a chain is fully specified by its matrix of
  mean-squared inter-residue distances, and conformations are drawn from it
  exactly: the bead coordinates are jointly Gaussian, with a covariance fixed by
  that matrix. The AFRC uses this, as do the self-avoiding walk models (as a
  Gaussian approximation with the right mean-squared distances).
* ``sample_freely_jointed_chain()``, ``sample_freely_rotating_chain()`` and
  ``sample_worm_like_chain()`` build chains bond by bond from the model's own
  rules, so they are exact (the worm-like chain up to a controlled
  discretization). Their exact mean-squared inter-residue distances are given by
  ``freely_rotating_chain_msd()`` and ``discrete_worm_like_chain_msd()``.

Ensembles are written as a PDB/XTC pair: the PDB holds the topology (one CA bead
per residue) and the first conformation, and the XTC holds every conformation.
Writing the XTC needs mdtraj, which is an optional dependency. Alternatively,
``write_multimodel_pdb()`` writes the whole ensemble as a single multi-model PDB
file, with no extra dependencies.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""

from __future__ import annotations

import os
import textwrap
from collections.abc import Callable, Iterator

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

from .exceptions import AFRCException

# one-letter to three-letter residue names, used for the PDB residue names
THREE_LETTER_CODES: dict[str, str] = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR',
}

# a (small) negative eigenvalue of the coordinate covariance is rounding noise; a
# large one means the distances cannot come from any Gaussian chain in 3D
_NEGATIVE_EIGENVALUE_TOLERANCE: float = 1e-8

# fixed-width PDB fields limit how large a structure (and ensemble) we can write
_MAX_PDB_RESIDUES: int = 9999
_MAX_PDB_MODELS: int = 9999
_MAX_PDB_COORDINATE: float = 999.999


# .....................................................................................
#
def gaussian_chain_factor(mean_squared_distances: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Factor the coordinate covariance of a Gaussian chain.

    For a Gaussian chain whose conformations are centred on the origin, each
    Cartesian coordinate of the beads is jointly Gaussian with covariance

    .. math::

       C = -\\frac{1}{6} J D J, \\qquad J = I - \\frac{1}{n}\\mathbf{1}\\mathbf{1}^T

    where :math:`D_{ij} = \\langle r_{ij}^2 \\rangle` is the matrix of
    mean-squared inter-residue distances. This returns a matrix :math:`A` with
    :math:`C = A A^T`, so that ``A @ z`` for a standard normal vector ``z`` is one
    coordinate (x, y or z) of one conformation.

    Parameters
    ----------
    mean_squared_distances : np.ndarray
        Symmetric [n x n] matrix of mean-squared inter-residue distances (in
        Angstroms squared), with zeros on the diagonal.

    Returns
    -------
    np.ndarray
        The [n x n] factor :math:`A`.

    Raises
    ------
    AFRCException
        If the matrix is not square and symmetric, or if the distances cannot be
        realized by any Gaussian chain in three dimensions (the covariance has a
        genuinely negative eigenvalue).

    """

    D = np.asarray(mean_squared_distances, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise AFRCException(f'The mean-squared distance matrix must be square (got shape {D.shape})')
    if not np.allclose(D, D.T):
        raise AFRCException('The mean-squared distance matrix must be symmetric')

    n = D.shape[0]
    if n == 0:
        return np.zeros((0, 0))
    J = np.eye(n) - 1.0/n
    covariance = -J @ D @ J / 6.0

    eigenvalues, eigenvectors = np.linalg.eigh(covariance)

    # the centring means one eigenvalue is exactly zero (up to rounding), and any
    # other tiny negative values are rounding noise too. A genuinely negative
    # eigenvalue means no Gaussian chain in 3D has these distances
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
    if np.min(eigenvalues) < -_NEGATIVE_EIGENVALUE_TOLERANCE * scale:
        raise AFRCException('These mean-squared distances cannot be realized by a Gaussian chain in three dimensions')

    return eigenvectors * np.sqrt(np.clip(eigenvalues, 0.0, None))


# .....................................................................................
#
def sample_gaussian_chain(factor: NDArray[np.float64], n_conformations: int, rng: np.random.Generator) -> NDArray[np.float64]:
    """
    Draw conformations of a Gaussian chain from its covariance factor.

    Parameters
    ----------
    factor : np.ndarray
        The [n x n] covariance factor from ``gaussian_chain_factor()``.

    n_conformations : int
        Number of conformations to draw.

    rng : np.random.Generator
        Random number generator to draw from.

    Returns
    -------
    np.ndarray
        Array of shape [n_conformations x n x 3] with the bead coordinates. Each
        conformation is centred on the origin and randomly oriented.

    """

    z = rng.standard_normal((n_conformations, factor.shape[0], 3))

    # factor @ z broadcasts over conformations: [n x n] @ [s x n x 3] -> [s x n x 3].
    # This goes through BLAS, and is ~25x faster than the equivalent einsum
    conformations: NDArray[np.float64] = factor @ z

    # the covariance is centred, so each conformation's centroid is zero in exact
    # arithmetic. When the smallest non-zero eigenvalues are nearly degenerate with
    # the centring's zero (as for small scaling exponents), rounding can leave a
    # ~1e-8 A offset; removing it changes no distance
    if conformations.shape[1] > 0:
        conformations -= conformations.mean(axis=1, keepdims=True)
    return conformations


# .....................................................................................
#
def mean_squared_distance_map(n_res: int, msd_of_separation: Callable[[NDArray[np.float64]], NDArray[np.float64]]) -> NDArray[np.float64]:
    """
    Build an [n x n] mean-squared distance matrix from a function of separation.

    For the homogeneous chain models (everything except the AFRC) the
    mean-squared distance between two residues depends only on how far apart
    they are in sequence, :math:`k = |i - j|`.

    Parameters
    ----------
    n_res : int
        Number of residues.

    msd_of_separation : callable
        Function mapping an array of separations :math:`k \\geq 1` to the
        mean-squared distances (in Angstroms squared).

    Returns
    -------
    np.ndarray
        Symmetric [n_res x n_res] matrix with zeros on the diagonal.

    """
    index = np.arange(n_res)
    separation = np.abs(np.subtract.outer(index, index)).astype(float)
    D = np.zeros((n_res, n_res))
    off_diagonal = separation > 0
    D[off_diagonal] = msd_of_separation(separation[off_diagonal])
    return D


# .....................................................................................
#
def freely_rotating_chain_msd(k: NDArray[np.float64], bond_length: float, cos_angle: float) -> NDArray[np.float64]:
    """
    Exact mean-squared distance between beads k bonds apart on a freely rotating chain.

    .. math::

       \\langle r^2 \\rangle = k b^2 \\frac{1 + \\alpha}{1 - \\alpha}
                           - 2 b^2 \\frac{\\alpha (1 - \\alpha^k)}{(1 - \\alpha)^2}

    where :math:`\\alpha` is the cosine of the angle between successive bonds.
    With :math:`\\alpha = 0` this is :math:`k b^2`, which is also exact for the
    freely jointed chain.

    Parameters
    ----------
    k : np.ndarray
        Number of bonds between the two beads.

    bond_length : float
        Bond length :math:`b` (in Angstroms).

    cos_angle : float
        :math:`\\alpha`, strictly between -1 and 1.

    Returns
    -------
    np.ndarray
        Mean-squared distances (in Angstroms squared).

    """
    k = np.asarray(k, dtype=float)
    b2 = bond_length * bond_length
    a = cos_angle
    msd: NDArray[np.float64] = k * b2 * (1 + a) / (1 - a) - 2 * b2 * a * (1 - np.power(a, k)) / (1 - a)**2
    return msd


# .....................................................................................
#
def worm_like_chain_msd(contour_length: NDArray[np.float64], persistence_length: float) -> NDArray[np.float64]:
    """
    Exact mean-squared end-to-end distance of a continuous worm-like chain.

    .. math::

       \\langle R^2 \\rangle = 2 L_p L - 2 L_p^2 \\left(1 - e^{-L/L_p}\\right)

    Parameters
    ----------
    contour_length : np.ndarray
        Contour length(s) :math:`L` (in Angstroms).

    persistence_length : float
        Persistence length :math:`L_p` (in Angstroms).

    Returns
    -------
    np.ndarray
        Mean-squared distances (in Angstroms squared).

    """
    L = np.asarray(contour_length, dtype=float)
    lp = persistence_length
    msd: NDArray[np.float64] = 2 * lp * L - 2 * lp * lp * (-np.expm1(-L / lp))
    return msd


# .....................................................................................
#
def discrete_worm_like_chain_msd(k: NDArray[np.float64], segment_length: float, subdivisions: int,
                                 tangent_correlation: float) -> NDArray[np.float64]:
    """
    Exact mean-squared distance between beads k residues apart on the discretized
    worm-like chain generated by ``sample_worm_like_chain()``.

    Each residue is split into ``subdivisions`` straight sub-segments of length
    :math:`s`, and successive sub-segments have tangent correlation :math:`c`, so
    this is the freely rotating chain result with :math:`k m` bonds of length
    :math:`s` and :math:`\\alpha = c`.

    Parameters
    ----------
    k : np.ndarray
        Separation in residues.

    segment_length : float
        Contour length per residue (in Angstroms).

    subdivisions : int
        Sub-segments per residue, :math:`m`.

    tangent_correlation : float
        Mean cosine between successive sub-segments, :math:`c`.

    Returns
    -------
    np.ndarray
        Mean-squared distances (in Angstroms squared).

    """
    s = segment_length / subdivisions
    return freely_rotating_chain_msd(np.asarray(k, dtype=float) * subdivisions, s, tangent_correlation)


# .....................................................................................
#
def _tangent_correlation(rule: str, sub_segment: float, persistence_length: float) -> float:
    """
    Tangent correlation between successive sub-segments of a discretized worm-like chain.

    Parameters
    ----------
    rule : str
        ``'exponential'`` for the worm-like chain's own correlation,
        :math:`e^{-s/L_p}`, or ``'matched'`` for :math:`(2L_p - s)/(2L_p + s)`,
        which makes the discretized chain's mean-squared size grow at exactly
        the continuous chain's rate, :math:`2 L_p` per unit contour length.

    sub_segment : float
        Sub-segment length :math:`s` (in Angstroms).

    persistence_length : float
        Persistence length :math:`L_p` (in Angstroms).

    Returns
    -------
    float
        The correlation :math:`c` (may be <= 0 for the matched rule if
        :math:`s \\geq 2 L_p`, which is not a usable discretization).

    """
    if rule == 'exponential':
        return float(np.exp(-sub_segment / persistence_length))
    return float((2 * persistence_length - sub_segment) / (2 * persistence_length + sub_segment))


# .....................................................................................
#
def worm_like_chain_discretization(segment_length: float, persistence_length: float, n_res: int,
                                   tolerance: float = 2e-4, max_subdivisions: int = 4096) -> tuple[int, float]:
    """
    Choose how to discretize the worm-like chain.

    Each residue is split into straight sub-segments, and we need the smallest
    number for which the discretized chain's mean-squared distances match the
    continuous worm-like chain to within ``tolerance`` (relative) at every
    separation in the chain. Two choices of the correlation between successive
    sub-segments are tried:

    * the worm-like chain's own, :math:`e^{-s/L_p}`, which is best for stiff
      chains;
    * :math:`(2L_p - s)/(2L_p + s)`, which reproduces the continuous chain's
      long-range size exactly and leaves only a short-range error of
      :math:`s^2/2`, so needs far fewer sub-segments when :math:`L_p` is small
      compared with a residue (e.g. 221 rather than 807 per residue at
      :math:`L_p = 0.1` A).

    Whichever needs fewer sub-segments is used.

    Parameters
    ----------
    segment_length : float
        Contour length per residue (in Angstroms).

    persistence_length : float
        Persistence length (in Angstroms).

    n_res : int
        Number of residues.

    tolerance : float
        Largest acceptable relative difference from the continuous chain.
        Default is 2e-4.

    max_subdivisions : int
        Upper limit on the number of sub-segments. Default is 4096.

    Returns
    -------
    tuple of (int, float)
        Sub-segments per residue, and the correlation between successive
        sub-segments.

    """
    if n_res < 2:
        return 1, _tangent_correlation('exponential', segment_length, persistence_length)
    k = np.arange(1, n_res, dtype=float)
    continuous = worm_like_chain_msd(k * segment_length, persistence_length)

    def error(rule: str, m: int) -> float:
        c = _tangent_correlation(rule, segment_length / m, persistence_length)
        if not 0 < c < 1:
            return float('inf')
        discrete = discrete_worm_like_chain_msd(k, segment_length, m, c)
        return float(np.max(np.abs(discrete / continuous - 1)))

    def smallest(rule: str) -> int:
        # double until good enough, then bisect back down to the smallest m that is
        upper = 1
        while error(rule, upper) > tolerance and upper < max_subdivisions:
            upper *= 2
        upper = min(upper, max_subdivisions)
        lower = upper // 2
        while upper - lower > 1:
            middle = (lower + upper) // 2
            if error(rule, middle) <= tolerance:
                upper = middle
            else:
                lower = middle
        return upper

    # prefer the worm-like chain's own correlation unless the matched one is cheaper
    m_exponential, m_matched = smallest('exponential'), smallest('matched')
    rule, m = ('matched', m_matched) if m_matched < m_exponential else ('exponential', m_exponential)
    return m, _tangent_correlation(rule, segment_length / m, persistence_length)


# .....................................................................................
#
def _unit_vectors(shape: tuple[int, ...], rng: np.random.Generator) -> NDArray[np.float64]:
    """
    Draw unit vectors uniformly on the sphere.

    Parameters
    ----------
    shape : tuple of int
        Leading shape; the result has shape ``shape + (3,)``.

    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Unit vectors.

    """
    v = rng.standard_normal(shape + (3,))
    unit: NDArray[np.float64] = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return unit


# .....................................................................................
#
def _rotate_away(t: NDArray[np.float64], cos_theta: NDArray[np.float64], phi: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Rotate unit vectors by a polar angle theta about a random azimuth phi.

    Returns :math:`\\cos\\theta\\, t + \\sin\\theta (\\cos\\phi\\, e_1 + \\sin\\phi\\, e_2)`,
    where :math:`e_1, e_2` complete an orthonormal basis with :math:`t`.

    Parameters
    ----------
    t : np.ndarray
        Unit vectors, shape [n x 3].

    cos_theta : np.ndarray
        Cosine of the angle between each input and output vector, shape [n].

    phi : np.ndarray
        Azimuthal angles, shape [n].

    Returns
    -------
    np.ndarray
        The rotated unit vectors, shape [n x 3].

    """
    # a reference axis that is never close to parallel with t
    reference = np.zeros_like(t)
    use_x = np.abs(t[:, 0]) < 0.9
    reference[use_x, 0] = 1.0
    reference[~use_x, 1] = 1.0

    e1 = np.cross(t, reference)
    e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(t, e1)

    sin_theta = np.sqrt(np.clip(1.0 - cos_theta * cos_theta, 0.0, None))
    out = (cos_theta[:, None] * t
           + (sin_theta * np.cos(phi))[:, None] * e1
           + (sin_theta * np.sin(phi))[:, None] * e2)

    # renormalize so rounding errors cannot accumulate along the chain
    rotated: NDArray[np.float64] = out / np.linalg.norm(out, axis=-1, keepdims=True)
    return rotated


# .....................................................................................
#
def _centre_chains(bonds: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Turn bond vectors into centred bead coordinates.

    Parameters
    ----------
    bonds : np.ndarray
        Bond vectors, shape [n_conformations x (n - 1) x 3].

    Returns
    -------
    np.ndarray
        Bead coordinates, shape [n_conformations x n x 3], each conformation
        centred on the origin.

    """
    positions = np.concatenate([np.zeros((bonds.shape[0], 1, 3)), np.cumsum(bonds, axis=1)], axis=1)
    centred: NDArray[np.float64] = positions - positions.mean(axis=1, keepdims=True)
    return centred


# .....................................................................................
#
def sample_freely_jointed_chain(n_res: int, bond_length: float, n_conformations: int, rng: np.random.Generator) -> NDArray[np.float64]:
    """
    Draw freely jointed chains: bonds of fixed length pointing in independent,
    uniformly random directions.

    Parameters
    ----------
    n_res : int
        Number of beads.

    bond_length : float
        Bond length (in Angstroms).

    n_conformations : int
        Number of conformations.

    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Bead coordinates, shape [n_conformations x n_res x 3], each conformation
        centred on the origin.

    """
    if n_res < 1:
        return np.zeros((n_conformations, 0, 3))
    bonds = bond_length * _unit_vectors((n_conformations, n_res - 1), rng)
    return _centre_chains(bonds)


# .....................................................................................
#
def sample_freely_rotating_chain(n_res: int, bond_length: float, cos_angle: float, n_conformations: int,
                                 rng: np.random.Generator) -> NDArray[np.float64]:
    """
    Draw freely rotating chains: bonds of fixed length with a fixed angle between
    successive bonds and uniformly random torsions.

    Parameters
    ----------
    n_res : int
        Number of beads.

    bond_length : float
        Bond length (in Angstroms).

    cos_angle : float
        Cosine of the angle between successive bond vectors.

    n_conformations : int
        Number of conformations.

    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Bead coordinates, shape [n_conformations x n_res x 3], each conformation
        centred on the origin.

    """
    if n_res < 1:
        return np.zeros((n_conformations, 0, 3))
    n_bonds = n_res - 1
    bonds = np.empty((n_conformations, n_bonds, 3))
    if n_bonds > 0:
        t = _unit_vectors((n_conformations,), rng)
        bonds[:, 0] = t
        cos_theta = np.full(n_conformations, cos_angle)
        for bond in range(1, n_bonds):
            t = _rotate_away(t, cos_theta, rng.uniform(0.0, 2 * np.pi, n_conformations))
            bonds[:, bond] = t
    return _centre_chains(bond_length * bonds)


# .....................................................................................
#
def _mean_cosine_to_concentration(mean_cosine: float) -> float:
    """
    Find the von Mises-Fisher concentration with a given mean cosine.

    Solves :math:`\\coth\\kappa - 1/\\kappa = c` for :math:`\\kappa`.

    Parameters
    ----------
    mean_cosine : float
        The target mean cosine :math:`c`, strictly between 0 and 1.

    Returns
    -------
    float
        The concentration :math:`\\kappa`.

    """
    def langevin(kappa: float) -> float:
        if kappa < 1e-4:
            return kappa / 3.0
        return float(1.0 / np.tanh(kappa) - 1.0 / kappa)

    return float(brentq(lambda kappa: langevin(kappa) - mean_cosine, 1e-12, 1e12, xtol=1e-14, rtol=1e-14))


# .....................................................................................
#
def _von_mises_fisher_cosines(kappa: float, size: tuple[int, ...], rng: np.random.Generator) -> NDArray[np.float64]:
    """
    Draw cos(theta) for directions from a von Mises-Fisher distribution in 3D.

    Uses the exact inverse-CDF form
    :math:`\\cos\\theta = 1 + \\log(u + (1 - u) e^{-2\\kappa}) / \\kappa`.

    Parameters
    ----------
    kappa : float
        Concentration (> 0).

    size : tuple of int
        Shape of the output.

    rng : np.random.Generator
        Random number generator.

    Returns
    -------
    np.ndarray
        Values of cos(theta) in [-1, 1].

    """
    u = rng.uniform(0.0, 1.0, size)
    cosines: NDArray[np.float64] = 1.0 + np.log(u + (1.0 - u) * np.exp(-2.0 * kappa)) / kappa
    return np.clip(cosines, -1.0, 1.0)


# .....................................................................................
#
def sample_worm_like_chain(n_res: int, segment_length: float, n_conformations: int, rng: np.random.Generator,
                           subdivisions: int, tangent_correlation: float) -> NDArray[np.float64]:
    """
    Draw worm-like chains.

    Each residue is split into ``subdivisions`` straight sub-segments of length
    :math:`s`. The direction of each sub-segment is drawn from a von Mises-Fisher
    distribution about the previous one, with the concentration set so the mean
    cosine between successive sub-segments is ``tangent_correlation``, and only
    the bead at the end of each residue is kept. Use
    ``worm_like_chain_discretization()`` to choose the discretization; the exact
    mean-squared distances of the resulting chain are
    ``discrete_worm_like_chain_msd()``.

    Parameters
    ----------
    n_res : int
        Number of beads.

    segment_length : float
        Contour length per residue (in Angstroms).

    n_conformations : int
        Number of conformations.

    rng : np.random.Generator
        Random number generator.

    subdivisions : int
        Sub-segments per residue.

    tangent_correlation : float
        Mean cosine between successive sub-segments, strictly between 0 and 1.

    Returns
    -------
    np.ndarray
        Bead coordinates, shape [n_conformations x n_res x 3], each conformation
        centred on the origin.

    Raises
    ------
    AFRCException
        If ``tangent_correlation`` is not strictly between 0 and 1.

    """
    if not 0 < tangent_correlation < 1:
        raise AFRCException(f'The tangent correlation must be strictly between 0 and 1 (got {tangent_correlation})')
    if n_res < 1:
        return np.zeros((n_conformations, 0, 3))

    s = segment_length / subdivisions
    kappa = _mean_cosine_to_concentration(tangent_correlation)

    # the tangents, held as separate x, y, z arrays so the inner loop is plain
    # elementwise arithmetic (np.cross and friends dominate the run time otherwise)
    t = _unit_vectors((n_conformations,), rng)
    tx, ty, tz = t[:, 0].copy(), t[:, 1].copy(), t[:, 2].copy()

    # random numbers are drawn a block of sub-segments at a time rather than per
    # sub-segment, with the block size capped so memory stays bounded for very
    # large ensembles
    block = max(1, min(subdivisions, (1 << 20) // max(n_conformations, 1)))

    n_bonds = n_res - 1
    bonds = np.zeros((n_conformations, n_bonds, 3))
    for bond in range(n_bonds):
        sx, sy, sz = np.zeros(n_conformations), np.zeros(n_conformations), np.zeros(n_conformations)
        for first in range(0, subdivisions, block):
            steps = min(block, subdivisions - first)
            cos_theta = _von_mises_fisher_cosines(kappa, (steps, n_conformations), rng)
            sin_theta = np.sqrt(np.clip(1.0 - cos_theta * cos_theta, 0.0, None))
            phi = rng.uniform(0.0, 2 * np.pi, (steps, n_conformations))
            along_e1, along_e2 = sin_theta * np.cos(phi), sin_theta * np.sin(phi)

            for step in range(steps):
                # the very first sub-segment keeps its random starting direction
                if bond > 0 or first + step > 0:
                    # e1 = normalize(t x ref), with ref the x axis unless t is close to it
                    # (t x x = (0, tz, -ty); t x y = (-tz, 0, tx)); e2 = t x e1
                    use_x = np.abs(tx) < 0.9
                    e1x = np.where(use_x, 0.0, -tz)
                    e1y = np.where(use_x, tz, 0.0)
                    e1z = np.where(use_x, -ty, tx)
                    inverse = 1.0 / np.sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
                    e1x, e1y, e1z = e1x * inverse, e1y * inverse, e1z * inverse
                    e2x, e2y, e2z = ty * e1z - tz * e1y, tz * e1x - tx * e1z, tx * e1y - ty * e1x

                    c, a, d = cos_theta[step], along_e1[step], along_e2[step]
                    nx, ny, nz = c * tx + a * e1x + d * e2x, c * ty + a * e1y + d * e2y, c * tz + a * e1z + d * e2z

                    # renormalize so rounding errors cannot accumulate along the chain
                    inverse = 1.0 / np.sqrt(nx * nx + ny * ny + nz * nz)
                    tx, ty, tz = nx * inverse, ny * inverse, nz * inverse
                sx += tx
                sy += ty
                sz += tz
        bonds[:, bond, 0], bonds[:, bond, 1], bonds[:, bond, 2] = s * sx, s * sy, s * sz
    return _centre_chains(bonds)


# .....................................................................................
#
def validate_n_conformations(n: object) -> int:
    """
    Check a requested number of conformations is a positive integer.

    Parameters
    ----------
    n : object
        The requested number (anything integer-valued, e.g. an ``int`` or
        ``np.int64``).

    Returns
    -------
    int
        The number, as an int.

    Raises
    ------
    AFRCException
        If ``n`` is not a positive integer.

    """
    try:
        count: int = int(n)  # type: ignore[call-overload]
    except (TypeError, ValueError):
        raise AFRCException(f'The number of conformations must be a positive integer (got {n})')
    if count != n or count < 1:
        raise AFRCException(f'The number of conformations must be a positive integer (got {n})')
    return count


# .....................................................................................
#
def validate_sequence(sequence: str) -> str:
    """
    Check a sequence can be written as a one-bead-per-residue structure.

    Parameters
    ----------
    sequence : str
        Amino acid sequence (case insensitive).

    Returns
    -------
    str
        The upper-case sequence.

    Raises
    ------
    AFRCException
        If the sequence is not a string, is empty, or contains a non-standard
        amino acid.

    """
    try:
        upper = sequence.upper()
    except AttributeError:
        raise AFRCException('The sequence must be a string of amino acids')
    if len(upper) == 0:
        raise AFRCException('The sequence must contain at least one amino acid')
    bad = sorted(set(upper) - set(THREE_LETTER_CODES))
    if bad:
        raise AFRCException(f'Sequence contains non-standard amino acids {bad}')
    return upper


# .....................................................................................
#
def save_conformations(conformations: NDArray[np.float64], sequence: str, filename: str, pdb_only: bool = False, remark: str = '') -> None:
    """
    Write conformations to disk as a PDB/XTC pair or a single multi-model PDB.

    This is the shared back end of every model's ``save_ensemble()``.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x n x 3] (in Angstroms).

    sequence : str
        Amino acid sequence of length n (case insensitive).

    filename : str
        Output path without an extension; ``.pdb`` (and ``.xtc``) are added. A
        trailing ``.pdb`` or ``.xtc`` is dropped first.

    pdb_only : bool
        If True, write one multi-model PDB instead of a PDB/XTC pair.

    remark : str
        Text for the ``REMARK`` lines of the PDB file.

    Raises
    ------
    AFRCException
        If the inputs are inconsistent or the files cannot be written (see
        ``write_ensemble()`` and ``write_multimodel_pdb()``).

    """
    sequence = validate_sequence(sequence)
    root, extension = os.path.splitext(filename)
    if extension.lower() in ('.pdb', '.xtc'):
        filename = root

    if pdb_only:
        write_multimodel_pdb(conformations, sequence, f'{filename}.pdb', remark=remark)
    else:
        write_ensemble(conformations, sequence, f'{filename}.pdb', f'{filename}.xtc',
                       remark=(remark + ', frame 1 shown') if remark else 'frame 1 shown')


# .....................................................................................
#
def _check_conformations(conformations: NDArray[np.float64], sequence: str) -> NDArray[np.float64]:
    """
    Check a set of conformations is consistent with a sequence.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x n x 3], or a single [n x 3] conformation.

    sequence : str
        One-letter amino acid sequence of length n.

    Returns
    -------
    np.ndarray
        The conformations as a float array of shape [n_conformations x n x 3].

    Raises
    ------
    AFRCException
        If the shapes do not match the sequence, the sequence is empty or too
        long for the PDB format, or it contains a non-standard amino acid.

    """

    xyz = np.asarray(conformations, dtype=float)
    if xyz.ndim == 2:
        xyz = xyz[np.newaxis]

    if len(sequence) == 0:
        raise AFRCException('Cannot write an ensemble for an empty sequence')
    if len(sequence) > _MAX_PDB_RESIDUES:
        raise AFRCException(f'The PDB format allows at most {_MAX_PDB_RESIDUES} residues (sequence has {len(sequence)})')
    if xyz.ndim != 3 or xyz.shape[1] != len(sequence) or xyz.shape[2] != 3 or xyz.shape[0] == 0:
        raise AFRCException(f'Conformations must have shape [n_conformations x {len(sequence)} x 3] (got {np.shape(conformations)})')

    bad = sorted(set(sequence) - set(THREE_LETTER_CODES))
    if bad:
        raise AFRCException(f'Sequence contains non-standard amino acids {bad}')

    return xyz


# .....................................................................................
#
def write_pdb(coordinates: NDArray[np.float64], sequence: str, filename: str, remark: str = '') -> None:
    """
    Write a single one-bead-per-residue conformation as a PDB file.

    Each residue is written as a single ``CA`` atom, with the residue name taken
    from the sequence, in chain A. Consecutive beads are joined with ``CONECT``
    records so that viewers draw the chain even though the beads are not at
    standard CA-CA spacing.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape [n x 3] with the bead coordinates, in Angstroms.

    sequence : str
        One-letter amino acid sequence of length n (upper-case).

    filename : str
        Path of the PDB file to write.

    remark : str
        Optional text written as a ``REMARK`` line at the top of the file.

    Raises
    ------
    AFRCException
        If the coordinates and sequence are inconsistent, or a coordinate is too
        large to fit in the fixed-width PDB format.

    """

    xyz = _check_conformations(coordinates, sequence)
    if xyz.shape[0] != 1:
        raise AFRCException('write_pdb() writes a single conformation; use write_multimodel_pdb() or write_ensemble() for several')

    _check_pdb_coordinates(xyz)

    lines = _remark_lines(remark) + _atom_lines(xyz[0], sequence) + _conect_lines(sequence) + ['END']
    with open(filename, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


# .....................................................................................
#
def write_multimodel_pdb(conformations: NDArray[np.float64], sequence: str, filename: str, remark: str = '') -> None:
    """
    Write a set of one-bead-per-residue conformations as a multi-model PDB file.

    Each conformation is written as its own ``MODEL``/``ENDMDL`` block, laid out
    as in ``write_pdb()``, so the whole ensemble lives in one file that mdtraj
    (and most viewers) read as a trajectory. This needs no extra dependencies,
    but the file is several times larger than a PDB/XTC pair, and the PDB format
    limits it to 9999 conformations.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x n x 3] with the bead coordinates, in
        Angstroms.

    sequence : str
        One-letter amino acid sequence of length n (upper-case).

    filename : str
        Path of the PDB file to write.

    remark : str
        Optional text written as ``REMARK`` lines at the top of the file.

    Raises
    ------
    AFRCException
        If the conformations and sequence are inconsistent, there are more than
        9999 conformations, or a coordinate is too large to fit in the
        fixed-width PDB format.

    """

    xyz = _check_conformations(conformations, sequence)
    if xyz.shape[0] > _MAX_PDB_MODELS:
        raise AFRCException(f'The PDB format allows at most {_MAX_PDB_MODELS} models (got {xyz.shape[0]}); write an XTC file instead')

    _check_pdb_coordinates(xyz)

    with open(filename, 'w') as fh:
        for chunk in _multimodel_pdb_chunks(xyz, sequence, remark):
            fh.write(chunk)


# .....................................................................................
#
def _multimodel_pdb_chunks(xyz: NDArray[np.float64], sequence: str, remark: str) -> Iterator[str]:
    """
    Generate the text of a multi-model PDB file, a piece at a time.

    The ``ATOM`` lines of every model are identical apart from their
    coordinates, so rather than formatting each line separately we build the
    block for one model once, with the coordinates left as ``%8.3f`` fields, and
    fill it in with a single ``%`` operation per model. This produces exactly the
    same text as formatting line by line (``write_pdb()``'s layout) but is ~4x
    faster, and yielding one model at a time means the whole file never has to
    be held in memory.

    Parameters
    ----------
    xyz : np.ndarray
        Validated coordinates, shape [n_conformations x n x 3] (Angstroms).

    sequence : str
        One-letter amino acid sequence of length n (upper-case).

    remark : str
        Text for the ``REMARK`` lines (may be empty).

    Yields
    ------
    str
        Consecutive pieces of the file; joined, they are the complete file.

    """
    for line in _remark_lines(remark):
        yield line + '\n'

    # the ATOM lines with their three coordinate columns (31-54) replaced by
    # format fields; the closing TER line has no coordinates. No other part of
    # these lines can contain a '%', so the template is safe to %-format
    atom_lines = _atom_lines(np.zeros((len(sequence), 3)), sequence)
    template = '\n'.join([line[:30] + '%8.3f%8.3f%8.3f' + line[54:] for line in atom_lines[:-1]] + [atom_lines[-1]])

    for model, frame in enumerate(xyz, start=1):
        yield f'MODEL     {model:4d}\n' + template % tuple(frame.ravel().tolist()) + '\nENDMDL\n'

    for line in _conect_lines(sequence):
        yield line + '\n'
    yield 'END\n'


# .....................................................................................
#
def _check_pdb_coordinates(xyz: NDArray[np.float64]) -> None:
    """
    Check every coordinate fits in the fixed-width PDB coordinate columns.

    Parameters
    ----------
    xyz : np.ndarray
        Coordinates (any shape), in Angstroms.

    Raises
    ------
    AFRCException
        If any coordinate's magnitude exceeds 999.999 A.

    """
    if np.max(np.abs(xyz)) > _MAX_PDB_COORDINATE:
        raise AFRCException(f'A coordinate exceeds {_MAX_PDB_COORDINATE} A, which does not fit in the PDB format; write an XTC file instead')


# .....................................................................................
#
def _remark_lines(remark: str) -> list[str]:
    """
    Format a remark as PDB ``REMARK`` lines.

    PDB lines are at most 80 characters, so a long remark is wrapped over
    several lines rather than cut off.

    Parameters
    ----------
    remark : str
        The remark text (may be empty).

    Returns
    -------
    list of str
        The ``REMARK`` lines (none for an empty remark).

    """
    return [f'REMARK   1 {chunk}' for chunk in textwrap.wrap(remark, width=69)]


# .....................................................................................
#
def _atom_lines(xyz: NDArray[np.float64], sequence: str) -> list[str]:
    """
    Format one conformation as PDB ``ATOM`` lines plus a closing ``TER``.

    Parameters
    ----------
    xyz : np.ndarray
        Array of shape [n x 3] with the bead coordinates, in Angstroms.

    sequence : str
        One-letter amino acid sequence of length n (upper-case).

    Returns
    -------
    list of str
        One ``ATOM`` line per residue (a ``CA`` bead in chain A), then ``TER``.

    """
    lines = []
    for index, (residue, (x, y, z)) in enumerate(zip(sequence, xyz), start=1):
        lines.append(f'ATOM  {index:5d}  CA  {THREE_LETTER_CODES[residue]} A{index:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           C')
    lines.append(f'TER   {len(sequence) + 1:5d}      {THREE_LETTER_CODES[sequence[-1]]} A{len(sequence):4d}')
    return lines


# .....................................................................................
#
def _conect_lines(sequence: str) -> list[str]:
    """
    Format ``CONECT`` records joining each bead to the next.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.

    Returns
    -------
    list of str
        One ``CONECT`` line per consecutive pair of beads.

    """
    return [f'CONECT{index:5d}{index + 1:5d}' for index in range(1, len(sequence))]


# .....................................................................................
#
def write_xtc(conformations: NDArray[np.float64], filename: str) -> None:
    """
    Write a set of conformations as an XTC trajectory.

    XTC files store coordinates in nanometres at a precision of 0.001 nm
    (0.01 A), so coordinates are converted from Angstroms on the way out.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x n x 3] with the bead coordinates, in
        Angstroms.

    filename : str
        Path of the XTC file to write.

    Raises
    ------
    AFRCException
        If mdtraj (needed to write XTC files) is not installed, or the array has
        the wrong shape.

    """

    try:
        from mdtraj.formats import XTCTrajectoryFile  # type: ignore[import-untyped]
    except ImportError:
        raise AFRCException('Writing XTC files needs mdtraj; install it with "pip install mdtraj" (or "pip install afrc[ensemble]")')

    xyz = np.asarray(conformations, dtype=float)
    if xyz.ndim != 3 or xyz.shape[2] != 3 or xyz.shape[0] == 0:
        raise AFRCException(f'Conformations must have shape [n_conformations x n x 3] (got {xyz.shape})')

    n_frames = xyz.shape[0]
    with XTCTrajectoryFile(filename, 'w') as fh:
        fh.write(xyz=(xyz/10.0).astype(np.float32),
                 time=np.arange(n_frames, dtype=np.float32),
                 step=np.arange(n_frames, dtype=np.int32))


# .....................................................................................
#
def write_ensemble(conformations: NDArray[np.float64], sequence: str, pdb_filename: str, xtc_filename: str, remark: str = '') -> None:
    """
    Write an ensemble as a PDB/XTC pair.

    The PDB file holds the topology and the first conformation; the XTC file
    holds every conformation (including the first). Load them together with,
    for example, ``mdtraj.load(xtc_filename, top=pdb_filename)`` or SOURSOP's
    ``SSTrajectory(xtc_filename, pdb_filename)``.

    Parameters
    ----------
    conformations : np.ndarray
        Array of shape [n_conformations x n x 3] with the bead coordinates, in
        Angstroms.

    sequence : str
        One-letter amino acid sequence of length n (upper-case).

    pdb_filename : str
        Path of the PDB file to write.

    xtc_filename : str
        Path of the XTC file to write.

    remark : str
        Optional text written as a ``REMARK`` line at the top of the PDB file.

    Raises
    ------
    AFRCException
        If the conformations and sequence are inconsistent, a coordinate in the
        first conformation is too large for the PDB format, or mdtraj is not
        installed.

    """

    xyz = _check_conformations(conformations, sequence)

    # check everything that could fail before writing either file, so a failure
    # never leaves half an ensemble behind. write_xtc() checks for mdtraj before
    # it creates its file, so writing the XTC first covers that case too
    _check_pdb_coordinates(xyz[0])

    write_xtc(xyz, xtc_filename)
    write_pdb(xyz[0], sequence, pdb_filename, remark=remark)
