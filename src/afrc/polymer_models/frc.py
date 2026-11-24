"""
frc.py

Freely rotating chain model with a tunable characteristic ratio.

Note that FRC here means freely *rotating* chain, not the Flory Random Coil
that the AFRC itself is built on.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from numpy.typing import NDArray
from afrc.config import P_OF_R_RESOLUTION
from afrc.ensemble import (freely_rotating_chain_msd, mean_squared_distance_map, sample_freely_rotating_chain,
                           save_conformations, validate_n_conformations)
from afrc.ensemble_report import EnsembleReport, ModelExpectations, compare_ensemble_to_model
from numpy.random import choice

class FRCException(Exception):
    """Exception raised by the freely rotating chain model."""
    pass

class FreelyRotatingChain:
    """
    Freely rotating chain of ``N`` bonds of length ``b`` with a fixed bond angle
    and free torsions.

    The freely rotating chain is an ideal chain: like the AFRC it has Gaussian
    end-to-end statistics with a scaling exponent of 0.5, but its size is set by
    a single stiffness parameter, the characteristic ratio :math:`C_\\infty`. The
    mean-squared end-to-end distance is the exact finite-N result

    .. math::

       \\langle R^2 \\rangle = C_\\infty N b^2 - 2 b^2 \\frac{\\alpha (1 - \\alpha^N)}{(1 - \\alpha)^2},
       \\qquad \\alpha = \\frac{C_\\infty - 1}{C_\\infty + 1}

    where :math:`\\alpha` is the cosine of the angle between successive bonds.
    ``c_inf = 1`` recovers the freely jointed chain, and ``c_inf = 2`` is a
    tetrahedral backbone. A freely rotating chain cannot reach the much larger
    characteristic ratio of a real polypeptide (:math:`C_\\infty \\approx 9`),
    which comes from hindered rotation - use the AFRC for that.

    This is a composition-independent reference model: the sequence is only used
    to set the number of bonds.

    References
    ----------
    [1] Flory, P. J. (1969). Statistical Mechanics of Chain Molecules.
    Wiley-Interscience.

    [2] Rubinstein, M., & Colby, R. H. (2003). Polymer Physics. Oxford
    University Press.

    """

    # .....................................................................................
    #
    def __init__(self, seq: str, p_of_r_resolution: float = P_OF_R_RESOLUTION, b: float = 3.8, c_inf: float = 2.0) -> None:
        """
        Create a FreelyRotatingChain object.

        Parameters
        ----------
        seq : str
            Amino acid sequence. Only its length (the number of bonds) is used,
            and it is not validated.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distribution. Default is
            0.05 A.

        b : float
            Bond length, in Angstroms. The default of 3.8 A is the Cα-Cα
            distance, i.e. one virtual bond per residue. Must be > 0.

        c_inf : float
            Characteristic ratio :math:`C_\\infty = (1 + \\alpha)/(1 - \\alpha)`,
            where :math:`\\alpha` is the cosine of the angle between successive
            bonds. ``c_inf = 1`` recovers the freely jointed chain; ``c_inf = 2``
            (the default) is a tetrahedral bond angle. Must be > 0.

        Raises
        ------
        FRCException
            If ``b`` or ``c_inf`` is not positive.

        """

        # cast and sanity check the free parameters
        self.b = float(b)
        if self.b <= 0:
            raise FRCException('Error, b (bond length) cannot be less than or equal to 0')

        self.c_inf = float(c_inf)
        if self.c_inf <= 0:
            raise FRCException('Error, c_inf (characteristic ratio) cannot be less than or equal to 0')

        # set sequence info - the number of bonds. The sequence itself is only used
        # to name residues when an ensemble is written to disk
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
            sum to 1).

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
        distribution. For long chains it approaches
        :math:`\\sqrt{C_\\infty} b \\sqrt{N}`.

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
        Build and cache the Gaussian end-to-end distribution.

        The mean-squared end-to-end distance is the exact finite-N freely
        rotating chain result (see the class docstring), and the distribution is
        the corresponding Gaussian

        .. math::

           P(r) = 4\\pi r^2 \\left( \\frac{3}{2\\pi \\langle R^2 \\rangle} \\right)^{3/2}
                  \\exp\\left( -\\frac{3 r^2}{2 \\langle R^2 \\rangle} \\right)

        evaluated from 0 to four times the root-mean-square size.

        """

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        N = self.nres
        b = self.b

        # the cosine of the angle between successive bonds, recovered from the
        # characteristic ratio: C∞ = (1 + α)/(1 - α)
        alpha = (self.c_inf - 1.0) / (self.c_inf + 1.0)

        # exact finite-N mean-squared end-to-end distance for the freely rotating chain.
        # The first term is the long-chain limit (C∞ N b²); the second is the finite-size
        # correction (which vanishes when α = 0, i.e. the freely jointed chain).
        mean_sq_re = self.c_inf * N * np.power(b, 2)
        if alpha != 0:
            mean_sq_re = mean_sq_re - 2 * np.power(b, 2) * alpha * (1 - np.power(alpha, N)) / np.power(1 - alpha, 2)

        # build an r-grid that comfortably captures the Gaussian tail regardless of the
        # chosen stiffness/segment length
        rms = np.sqrt(mean_sq_re)
        p_dist = np.arange(0, 4.0*rms, self.p_of_r_resolution)

        # standard Gaussian chain end-to-end distribution
        A = np.power(3.0/(2*np.pi*mean_sq_re), 1.5)
        p_val_raw = 4*np.pi*np.power(p_dist, 2)*A*np.exp(-(3*np.power(p_dist, 2))/(2*mean_sq_re))

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw / np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


    # .....................................................................................
    #
    def _cos_angle(self) -> float:
        """
        Return the cosine of the angle between successive bond vectors.

        Returns
        -------
        float
            :math:`\\alpha = (C_\\infty - 1)/(C_\\infty + 1)`.

        """
        return (self.c_inf - 1.0) / (self.c_inf + 1.0)


    # .....................................................................................
    #
    def get_mean_squared_distance_map(self) -> NDArray[np.float64]:
        """
        Return the exact mean-squared distance between every pair of residues.

        For beads k bonds apart this is the exact finite-k freely rotating chain
        result (see the class docstring) - the same expression the model uses for
        the whole chain.

        Returns
        -------
        np.ndarray
            Symmetric [N x N] matrix of mean-squared distances (in Angstroms
            squared), with zeros on the diagonal.

        """
        cos_angle = self._cos_angle()
        return mean_squared_distance_map(self.nres, lambda k: freely_rotating_chain_msd(k, self.b, cos_angle))


    # .....................................................................................
    #
    def sample_conformations(self, n: int = 1000, seed: int | np.random.Generator | None = None) -> NDArray[np.float64]:
        """
        Generate 3D conformations (one bead per residue) of the freely rotating chain.

        Consecutive beads are joined by bonds of exactly ``b``, every pair of
        consecutive bonds meets at the same bond angle (set by ``c_inf``), and
        each torsion is uniformly random - the freely rotating chain itself. Its
        mean-squared distances are exact; its distance distributions are not
        Gaussian for short separations, unlike the model's analytical P(r),
        which is a Gaussian with the right mean-squared size. Note that
        ``c_inf = 1`` gives a fixed 90 degree bond angle, which has the same
        mean-squared size as the freely jointed chain but is a different chain.

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
        return sample_freely_rotating_chain(self.nres, self.b, self._cos_angle(), validate_n_conformations(n), np.random.default_rng(seed))


    # .....................................................................................
    #
    def save_ensemble(self, filename: str, n: int = 1000, seed: int | np.random.Generator | None = None,
                      pdb_only: bool = False) -> NDArray[np.float64]:
        """
        Generate a freely rotating chain ensemble and write it to disk.

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
                           remark=(f'Freely rotating chain ensemble (b = {self.b:g} A, c_inf = {self.c_inf:g}) from afrc {__version__}: '
                                   f'{len(conformations)} conformations'))
        return conformations


    # .....................................................................................
    #
    def check_ensemble(self, conformations: NDArray[np.float64]) -> EnsembleReport:
        """
        Check how well an ensemble reproduces the freely rotating chain.

        The exact checks are the mean-squared distance of every residue pair (in
        bands of separation), the root-mean-square radius of gyration and
        first-to-last distance that follow from them, and the chain's geometry:
        every bond exactly ``b`` long, every bond angle fixed by ``c_inf``, and
        no pair further apart than its contour length. The model's analytical
        (Gaussian) end-to-end distribution and its radius of gyration are
        reported as context.

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
            end_to_end = FreelyRotatingChain('A' * (self.nres - 1), self.p_of_r_resolution, b=self.b,
                                             c_inf=self.c_inf).get_end_to_end_distribution()

        expectations = ModelExpectations(
            'Freely rotating chain', self.get_mean_squared_distance_map(),
            bond_length=self.b, bond_angle_cosine=self._cos_angle(), contour_length_per_residue=self.b,
            reference_rg=self.get_mean_radius_of_gyration(),
            reference_rg_label='RMS Rg (sqrt(<Re^2>/6) for N bonds)',
            reference_re=self.get_root_mean_squared_end_to_end_distance(),
            reference_re_label='whole-chain RMS Re (N bonds)',
            end_to_end_distribution=end_to_end,
            end_to_end_label='the Gaussian P(r) with the exact <R^2> for N-1 bonds (Gaussian only for long chains)')
        return compare_ensemble_to_model(conformations, expectations)

