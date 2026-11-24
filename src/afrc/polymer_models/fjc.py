"""
fjc.py

Freely jointed chain (FJC) model using the non-Gaussian Kuhn-Grün distribution.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from numpy.typing import NDArray
from afrc.config import P_OF_R_RESOLUTION
from afrc.ensemble import (freely_rotating_chain_msd, mean_squared_distance_map, sample_freely_jointed_chain,
                           save_conformations, validate_n_conformations)
from afrc.ensemble_report import EnsembleReport, ModelExpectations, compare_ensemble_to_model
from numpy.random import choice

class FJCException(Exception):
    """Exception raised by the freely jointed chain model."""
    pass

class FreelyJointedChain:
    """
    Freely jointed chain of ``N`` rigid segments of length ``b``.

    Unlike the (Gaussian) AFRC, the end-to-end distribution here is the
    non-Gaussian Kuhn-Grün distribution [1], which respects the finite
    extensibility of the chain - the end-to-end distance can never exceed the
    contour length :math:`L = Nb`. At small extensions it reduces to the Gaussian
    result, so for typical IDR lengths the bulk of the distribution is close to
    the AFRC and the differences appear in the tail.

    This is a composition-independent reference model: the sequence is only used
    to set the number of segments.

    References
    ----------
    [1] Kuhn, W., & Grün, F. (1942). Beziehungen zwischen elastischen Konstanten
    und Dehnungsdoppelbrechung hochelastischer Stoffe. Kolloid-Zeitschrift,
    101(3), 248-271.

    [2] Cohen, A. (1991). A Padé approximant to the inverse Langevin function.
    Rheologica Acta, 30(3), 270-273.

    """

    # .....................................................................................
    #
    def __init__(self, seq: str, p_of_r_resolution: float = P_OF_R_RESOLUTION, b: float = 3.8) -> None:
        """
        Create a FreelyJointedChain object.

        Parameters
        ----------
        seq : str
            Amino acid sequence. Only its length (the number of segments) is
            used, and it is not validated.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distribution. Default is
            0.05 A.

        b : float
            Segment (Kuhn) length, in Angstroms. The default of 3.8 A is the Cα-Cα
            distance, i.e. one segment per residue. Must be > 0.

        Raises
        ------
        FJCException
            If ``b`` is not positive.

        """

        # cast and sanity check
        self.b = float(b)
        if self.b <= 0:
            raise FJCException('Error, b (segment length) cannot be less than or equal to 0')

        # set sequence info - the number of segments. The sequence itself is only
        # used to name residues when an ensemble is written to disk
        self.nres = len(seq)
        self.seq = seq

        # p_of_r_resolution defines the P(r) resolution in angstroms - i.e. basically
        # the spacing between r values in a P(r) vs. r plot
        self.p_of_r_resolution = p_of_r_resolution

        # set distribution info to false - these are calculated if/when needed
        self.__p_of_Re_R = False
        self.__p_of_Re_P = False

        # this sets a flag that is useful for letting certain functions work when
        # there's a chain length of 0
        if len(seq) == 0:
            self.zero_length = True
        else:
            self.zero_length = False


    # .....................................................................................
    #
    def get_end_to_end_distribution(self):
        """
        Return the end-to-end distance distribution.

        The distribution is computed on first use and then cached.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1). The grid never extends past the contour length.

        """
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self):
        """
        Return the mean end-to-end distance, :math:`\\langle R_e \\rangle`.

        This is the expectation over the end-to-end distribution,
        :math:`\\sum r P(r)`.

        Returns
        -------
        float
            The mean end-to-end distance (in Angstroms).

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sum(a * b)


    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Return the root-mean-square end-to-end distance,
        :math:`\\sqrt{\\langle R_e^2 \\rangle}`.

        This is the square root of :math:`\\sum r^2 P(r)` over the end-to-end
        distribution. For long chains it approaches the ideal-chain value
        :math:`b\\sqrt{N}` from below.

        Returns
        -------
        float
            The root-mean-square end-to-end distance (in Angstroms).

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sqrt(np.sum(b * np.power(a, 2)))


    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self):
        """
        Return the root-mean-square radius of gyration,
        :math:`\\sqrt{\\langle R_g^2 \\rangle}`.

        This uses the ideal-chain relation
        :math:`\\langle R_g^2 \\rangle = \\langle R_e^2 \\rangle / 6`. Note that
        despite the method name this is the *root-mean-square* radius of
        gyration, not :math:`\\langle R_g \\rangle`. The name is kept for
        consistency with the other models.

        Returns
        -------
        float
            The root-mean-square radius of gyration (in Angstroms).

        """
        return self.get_root_mean_squared_end_to_end_distance() / np.sqrt(6)


    # .....................................................................................
    #
    def sample_end_to_end_distribution(self, n=1000):
        """
        Draw random end-to-end distances from the distribution.

        Useful for building a size-matched sample to compare against simulation
        data.

        Parameters
        ----------
        n : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``n`` independent end-to-end distances (in Angstroms). For a
            zero-length chain every value is 0.

        """
        if self.zero_length:
            return np.repeat(0.0, n)

        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return choice(self.__p_of_Re_R, n, p=self.__p_of_Re_P)


    # .....................................................................................
    #
    def __compute_end_to_end_distribution(self):
        """
        Build and cache the Kuhn-Grün end-to-end distribution.

        .. math::

           P(r) \\propto 4\\pi r^2 \\exp\\left[ -N \\left( x\\beta +
                  \\ln\\frac{\\beta}{\\sinh\\beta} \\right) \\right],
           \\qquad x = \\frac{r}{Nb}

        where :math:`\\beta = \\mathcal{L}^{-1}(x)` is the inverse Langevin
        function, evaluated with the Cohen Padé approximant
        :math:`\\beta \\approx x(3 - x^2)/(1 - x^2)`. The distribution is only
        defined for :math:`r < Nb`, so the grid stops short of the contour length.

        """

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        # contour length
        L = self.nres * self.b

        # use the same style of r-grid as the other models, but make sure it always
        # reaches four times the ideal-chain size b*sqrt(N) - the fixed 21*sqrt(N)
        # grid is fine for b = 3.8 A but cut into the tail for larger segment lengths
        # (by b = 20 A it truncated the RMS by over 25%). Never exceed the contour
        # length, where the FJC distribution is undefined
        r_upper = min(max(3*(7*np.power(self.nres, 0.5)), 4*self.b*np.power(self.nres, 0.5)), L)
        p_dist = np.arange(0, r_upper, self.p_of_r_resolution)

        # fractional extension (strictly < 1 because arange excludes the endpoint)
        x = p_dist / L

        # the Cohen Padé approximant to the inverse Langevin function, and a
        # numerically stable evaluation of ln(β / sinh β) that does not overflow
        # for large β (i.e. as the chain approaches full extension)
        with np.errstate(divide='ignore', invalid='ignore'):
            beta = x * (3 - np.power(x, 2)) / (1 - np.power(x, 2))
            ln_sinh_beta = beta + np.log1p(-np.exp(-2*beta)) - np.log(2)
            ln_beta_over_sinh = np.log(beta) - ln_sinh_beta

            exponent = -self.nres * (x*beta + ln_beta_over_sinh)
            p_val_raw = np.power(p_dist, 2) * np.exp(exponent)

        # the r = 0 point evaluates to 0/0; the r^2 prefactor makes P(0) = 0 anyway
        p_val_raw[0] = 0.0
        p_val_raw = np.nan_to_num(p_val_raw, nan=0.0, posinf=0.0, neginf=0.0)

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw / np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


    # .....................................................................................
    #
    def get_mean_squared_distance_map(self) -> NDArray[np.float64]:
        """
        Return the exact mean-squared distance between every pair of residues.

        Beads k residues apart are joined by k freely jointed bonds, so
        :math:`\\langle r^2 \\rangle = k b^2` exactly.

        Returns
        -------
        np.ndarray
            Symmetric [N x N] matrix of mean-squared distances (in Angstroms
            squared), with zeros on the diagonal.

        """
        return mean_squared_distance_map(self.nres, lambda k: freely_rotating_chain_msd(k, self.b, 0.0))


    # .....................................................................................
    #
    def sample_conformations(self, n: int = 1000, seed: int | np.random.Generator | None = None) -> NDArray[np.float64]:
        """
        Generate 3D conformations (one bead per residue) of the freely jointed chain.

        Consecutive beads are joined by bonds of exactly ``b`` pointing in
        independent, uniformly random directions - the freely jointed chain
        itself, not an approximation to it. (Note that the model's analytical
        end-to-end distribution, the Kuhn-Grün form, is itself an approximation
        that becomes exact for long chains.)

        Parameters
        ----------
        n : int
            Number of conformations. Default is 1000.

        seed : int, np.random.Generator or None
            Seed (or generator) for reproducible ensembles. If None (default) a
            fresh, unpredictable seed is used.

        Returns
        -------
        np.ndarray
            Array of shape [n x N x 3] with the bead coordinates, in Angstroms,
            each conformation centred on the origin.

        Raises
        ------
        AFRCException
            If ``n`` is not a positive integer.

        """
        return sample_freely_jointed_chain(self.nres, self.b, validate_n_conformations(n), np.random.default_rng(seed))


    # .....................................................................................
    #
    def save_ensemble(self, filename: str, n: int = 1000, seed: int | np.random.Generator | None = None,
                      pdb_only: bool = False) -> NDArray[np.float64]:
        """
        Generate a freely jointed chain ensemble and write it to disk.

        Writes ``<filename>.pdb`` (topology and first conformation) and
        ``<filename>.xtc`` (every conformation; needs mdtraj), or with
        ``pdb_only=True`` a single multi-model ``<filename>.pdb``. See
        ``sample_conformations()`` for how the conformations are generated.

        Parameters
        ----------
        filename : str
            Output path without an extension (a trailing ``.pdb`` or ``.xtc`` is
            dropped).

        n : int
            Number of conformations. Default is 1000.

        seed : int, np.random.Generator or None
            Seed (or generator) for reproducible ensembles.

        pdb_only : bool
            Write a single multi-model PDB instead of a PDB/XTC pair. Default is
            False.

        Returns
        -------
        np.ndarray
            The conformations that were written, shape [n x N x 3] (Angstroms).

        Raises
        ------
        AFRCException
            If ``n`` is invalid, the sequence has non-standard amino acids, mdtraj
            is missing (for a PDB/XTC pair), or the PDB format limits are exceeded.

        """
        from afrc import __version__

        conformations = self.sample_conformations(n=n, seed=seed)
        save_conformations(conformations, self.seq, filename, pdb_only=pdb_only,
                           remark=f'Freely jointed chain ensemble (b = {self.b:g} A) from afrc {__version__}: {len(conformations)} conformations')
        return conformations


    # .....................................................................................
    #
    def check_ensemble(self, conformations: NDArray[np.float64]) -> EnsembleReport:
        """
        Check how well an ensemble reproduces the freely jointed chain.

        The exact checks are the mean-squared distance of every residue pair
        (:math:`k b^2`, in bands of separation), the root-mean-square radius of
        gyration and first-to-last distance that follow from them, and the
        chain's geometry: every bond exactly ``b`` long and no pair further apart
        than its contour length. The model's analytical (Kuhn-Grün) end-to-end
        distribution and its radius of gyration are reported as context.

        Parameters
        ----------
        conformations : np.ndarray
            Array of shape [n_conformations x N x 3] (in Angstroms).

        Returns
        -------
        EnsembleReport
            The checks; ``report.passed`` says whether the ensemble is model-like
            and ``report.format()`` gives a printable report.

        Raises
        ------
        AFRCException
            If the conformations do not match this chain.

        """
        end_to_end = None
        if self.nres >= 2:
            end_to_end = FreelyJointedChain('A' * (self.nres - 1), self.p_of_r_resolution, b=self.b).get_end_to_end_distribution()

        expectations = ModelExpectations(
            'Freely jointed chain', self.get_mean_squared_distance_map(),
            bond_length=self.b, contour_length_per_residue=self.b,
            reference_rg=self.get_mean_radius_of_gyration(),
            reference_rg_label='RMS Rg (sqrt(<Re^2>/6) for N segments)',
            reference_re=self.get_root_mean_squared_end_to_end_distance(),
            reference_re_label='whole-chain RMS Re (N segments)',
            end_to_end_distribution=end_to_end,
            end_to_end_label='the Kuhn-Grün P(r) for N-1 segments (an approximation that is exact for long chains)')
        return compare_ensemble_to_model(conformations, expectations)

