"""
wlc2.py

Worm-like chain (WLC) model using the closed form of O'Brien et al. (2009).

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from numpy.typing import NDArray
from afrc.ensemble import (discrete_worm_like_chain_msd, mean_squared_distance_map, sample_worm_like_chain, save_conformations,
                           validate_n_conformations, worm_like_chain_discretization, worm_like_chain_msd)
from afrc.ensemble_report import EnsembleReport, ModelExpectations, compare_ensemble_to_model
from afrc.config import P_OF_R_RESOLUTION

class WLC2Exception(Exception):
    """Exception raised by the O'Brien worm-like chain model."""
    pass

class WormLikeChain2:
    """
    Worm-like chain model, as implemented by O'Brien et al. (2009).

    This is a composition-independent reference model: the sequence is only used
    to set the number of residues, and hence the contour length
    :math:`L_c = N b`. Unlike the Zhou model (``WormLikeChain``) the O'Brien
    expression enforces finite extensibility exactly and stays well behaved for
    long chains, and this model also provides a closed-form radius of gyration.

    References
    ----------
    [1] O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009).
    How accurate are polymer models in the analysis of Forster resonance
    energy transfer experiments on proteins? The Journal of Chemical Physics,
    130(12), 124903.

    """

    # .....................................................................................
    #
    def __init__(self, seq: str, p_of_r_resolution: float = P_OF_R_RESOLUTION, lp: float = 3.0, aa_size: float = 3.8) -> None:
        """
        Create a WormLikeChain2 object.

        Parameters
        ----------
        seq : str
            Amino acid sequence. Only its length is used, and it is not
            validated.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distribution. Default is
            0.05 A.

        lp : float
            Persistence length, in Angstroms. We use a default of 3.0 A, although
            4 A is also common in the literature. Must be > 0.

        aa_size : float
            Contour length per residue (called :math:`b` in the literature), in
            Angstroms. The default of 3.8 A is the Cα-Cα distance. Must be > 0.

        Raises
        ------
        WLC2Exception
            If ``lp`` or ``aa_size`` is not positive, or if the contour length
            (``len(seq) * aa_size``) is shorter than the persistence length.

        """

        # set sequence info. The sequence itself is only used to name residues when
        # an ensemble is written to disk
        self.nres = len(seq)
        self.seq = seq
        self._discretization = None

        # first cast to floats
        self.lp = float(lp)
        self.b = float(aa_size)

        # also sanity check

        if self.lp <= 0:
            raise WLC2Exception('Error, lp cannot be less than or equal to 0')

        if self.b <= 0:
            raise WLC2Exception('Error, aa_size cannot be less than or equal to 0')

        Lc = self.b*self.nres

        # the chain must be at least one persistence length long. Note this compares
        # the CONTOUR length (N*aa_size, in Angstroms) with the persistence length -
        # this check previously compared the number of residues against lp, which mixes
        # a residue count with a length in Angstroms
        if Lc < self.lp:
            raise WLC2Exception('Passed sequence has a contour length (%.2f A) shorter than the persistence length (%.2f A)' % (Lc, self.lp))

        # next calculate params as defined by O'Brien et al
        self.alpha = (3*Lc) / (4*self.lp)
        self.C2 = 1/(2*self.lp)

        # C1 is the analytical normalization constant. Note that exp(alpha) overflows
        # for alpha > ~709 (about 750 residues at the default lp), so we build it in
        # log space and let it go to inf quietly for very long chains. C1 is kept for
        # reference only - the distribution itself is evaluated in log space and
        # normalized numerically, so it never depends on this value
        log_C1 = self.alpha + 1.5*np.log(self.alpha) - 1.5*np.log(np.pi) - np.log(1 + 3/self.alpha + (15/4)/np.power(self.alpha, 2))
        with np.errstate(over='ignore'):
            self.C1 = np.exp(log_C1)

        # p_of_r_resolution defines the P(r) resolution in angstroms - i.e. basically
        # the spacing between r values in a P(r) vs. r plot
        self.p_of_r_resolution = p_of_r_resolution

        # set distribution info to false - these are calculated if/when needed
        self.__p_of_Re_R = False
        self.__p_of_Re_P = False

        # an empty sequence has zero contour length and so has already failed the
        # contour-length check above; the flag is kept for interface parity with the
        # other models
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
            sum to 1). No probability lies at or beyond the contour length.

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
        [a,b] = self.get_end_to_end_distribution()

        return np.sum(a*b)


    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Return the root-mean-square end-to-end distance,
        :math:`\\sqrt{\\langle R_e^2 \\rangle}`.

        This is the square root of :math:`\\sum r^2 P(r)` over the end-to-end
        distribution.

        Returns
        -------
        float
            The root-mean-square end-to-end distance (in Angstroms).

        """


        [a,b] = self.get_end_to_end_distribution()
        return np.sqrt(np.sum(b*np.power(a,2)))



    # .....................................................................................
    #
    def __compute_end_to_end_distribution(self):
        """
        Build and cache the end-to-end distribution.

        With :math:`x = r/L_c` and :math:`\\alpha = 3L_c/(4L_p)`, the O'Brien
        expression is

        .. math::

           P(r) = \\frac{4\\pi C_1 r^2}{L_c^3 (1 - x^2)^{9/2}}
                  \\exp\\left( -\\frac{\\alpha}{1 - x^2} \\right)

        We evaluate its logarithm and normalize numerically, which avoids the
        overflow in :math:`C_1 \\propto e^{\\alpha}` for long chains. The grid
        runs from 0 to the smaller of :math:`L_c` and
        :math:`\\max(21\\sqrt{N}, 4\\sqrt{2 L_p L_c})`.

        """

        # define persistence length and contour length
        Lp = self.lp
        Lc = self.nres*self.b

        # use the same r-grid as the parent AFRC model, but make sure it always
        # reaches four times the ideal-chain size sqrt(2*Lp*Lc) - the previous grid
        # depended on Lp but not on aa_size, and cut into the tail for stiff chains
        # or large segment sizes. Nothing lies beyond the contour length, so the grid
        # never needs to go past it
        upper = min(max(3*(7*np.power(self.nres, 0.5)), 4*np.sqrt(2*Lp*Lc)), Lc)
        p_dist = np.arange(0, upper, self.p_of_r_resolution)

        # evaluate log P(r) (up to a constant) across the grid:
        #
        #   log P = 2 log r - (9/2) log(1 - x^2) - alpha/(1 - x^2),   x = r/Lc
        #
        # where alpha = 3*Lc/(4*Lp). Working in log space matters for long chains:
        # the linear-space form multiplies C1 ~ exp(alpha) by exp(-alpha/(1 - x^2)),
        # which overflows and underflows respectively once alpha passes ~709 and
        # previously returned an all-NaN distribution. At r = 0 the log is -inf,
        # which correctly gives P(0) = 0
        x2 = np.power(p_dist/Lc, 2)
        with np.errstate(divide='ignore', invalid='ignore'):
            log_p = 2*np.log(p_dist) - 4.5*np.log(1 - x2) - self.alpha/(1 - x2)

        # the grid stops short of Lc, but guard against a floating-point grid point
        # landing on (or past) it - the chain can never be longer than its contour
        # length, so such points carry zero probability
        log_p[x2 >= 1] = -np.inf

        # shift by the maximum before exponentiating so the largest value is 1
        p_val_raw = np.exp(log_p - np.max(log_p))

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self):
        """
        Return the root-mean-square radius of gyration,
        :math:`\\sqrt{\\langle R_g^2 \\rangle}`.

        O'Brien et al. give :math:`\\langle R_g^2 \\rangle` in closed form (they
        do not say so explicitly, but the expression is a mean-square). With
        :math:`C_2 = 1/(2L_p)`,

        .. math::

           \\langle R_g^2 \\rangle = \\frac{L_c}{6 C_2} - \\frac{1}{4 C_2^2}
                + \\frac{1}{4 C_2^3 L_c} - \\frac{1 - e^{-L_c/L_p}}{8 C_2^4 L_c^2}

        which is the Benoit-Doty worm-like chain result
        :math:`L_c L_p/3 - L_p^2 + 2L_p^3/L_c - 2L_p^4/L_c^2 (1 - e^{-L_c/L_p})`.
        It reduces to :math:`L_c L_p/3` for a flexible chain and to the rigid-rod
        value :math:`L_c^2/12` when :math:`L_c \\ll L_p`.

        Note that despite the method name this is the *root-mean-square* radius
        of gyration, not :math:`\\langle R_g \\rangle`. The name is kept for
        consistency with the other models.

        Returns
        -------
        float
            The root-mean-square radius of gyration (in Angstroms).

        References
        ----------
        [1] O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009).
        How accurate are polymer models in the analysis of Forster resonance
        energy transfer experiments on proteins? The Journal of Chemical Physics,
        130(12), 124903.

        [2] Benoit, H., & Doty, P. (1953). Light scattering from non-Gaussian
        chains. The Journal of Physical Chemistry, 57(9), 958-963.

        """

        Lc = self.nres*self.b
        C2 = self.C2
        Lp = self.lp

        # NOTE: the second term carries a MINUS sign (it is -Lp^2 in the
        # Benoit-Doty form); it was previously written as +1/(4*C2^2), which
        # over-estimated Rg (badly so for short chains, where it left a
        # spurious constant 2*Lp^2 offset instead of the correct Lc^2/12
        # rigid-rod limit)
        return np.sqrt(Lc/(6*C2) - 1/(4*np.power(C2,2)) +  1/(Lc*4*np.power(C2, 3)) - (1 - np.exp(-Lc/Lp))/(8*np.power(C2, 4)*np.power(Lc, 2)))


    # .....................................................................................
    #
    def _ensemble_discretization(self) -> tuple[int, float]:
        """
        Return how the chain is discretized to generate ensembles.

        Chosen (once, then cached) by ``worm_like_chain_discretization()`` so the
        discretized chain matches the continuous worm-like chain's mean-squared
        distances to within 0.02% at every separation in the chain.

        Returns
        -------
        tuple of (int, float)
            Sub-segments per residue, and the correlation between successive
            sub-segments.

        """
        if self._discretization is None:
            self._discretization = worm_like_chain_discretization(self.b, self.lp, self.nres)
        return self._discretization


    # .....................................................................................
    #
    def get_mean_squared_distance_map(self) -> NDArray[np.float64]:
        """
        Return the exact mean-squared distance between every pair of residues.

        For beads k residues apart (contour length :math:`L = k b`) this is the
        exact worm-like chain result
        :math:`\\langle r^2 \\rangle = 2 L_p L - 2 L_p^2 (1 - e^{-L/L_p})`.

        Returns
        -------
        np.ndarray
            Symmetric [N x N] matrix of mean-squared distances (in Angstroms
            squared), with zeros on the diagonal.

        """
        return mean_squared_distance_map(self.nres, lambda k: worm_like_chain_msd(k * self.b, self.lp))


    # .....................................................................................
    #
    def sample_conformations(self, n: int = 1000, seed: int | np.random.Generator | None = None) -> NDArray[np.float64]:
        """
        Generate 3D conformations (one bead per residue) of the worm-like chain.

        Each residue is split into short straight sub-segments, each bending away
        from the previous one with a fixed mean cosine, and only the bead at the
        end of each residue is kept. That correlation is either the worm-like
        chain's own, :math:`e^{-s/L_p}` for sub-segments of length :math:`s`, or -
        for chains that are flexible on the scale of a residue - one that
        reproduces its long-range size exactly with far fewer sub-segments (see
        ``afrc.ensemble.worm_like_chain_discretization()``). Either way the
        number of sub-segments is chosen so the mean-squared distances match the
        continuous worm-like chain to within 0.02% at every separation. This
        samples the worm-like chain itself; the model's analytical end-to-end
        distribution (O'Brien 2009) is an approximation to it.

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
        m, correlation = self._ensemble_discretization()
        return sample_worm_like_chain(self.nres, self.b, validate_n_conformations(n), np.random.default_rng(seed), m, correlation)


    # .....................................................................................
    #
    def save_ensemble(self, filename: str, n: int = 1000, seed: int | np.random.Generator | None = None,
                      pdb_only: bool = False) -> NDArray[np.float64]:
        """
        Generate a worm-like chain ensemble and write it to disk.

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
                           remark=(f'Worm-like chain ensemble (lp = {self.lp:g} A, aa_size = {self.b:g} A) from afrc {__version__}: '
                                   f'{len(conformations)} conformations'))
        return conformations


    # .....................................................................................
    #
    def check_ensemble(self, conformations: NDArray[np.float64]) -> EnsembleReport:
        """
        Check how well an ensemble reproduces the worm-like chain.

        The exact checks are the mean-squared distance of every residue pair (in
        bands of separation), the root-mean-square radius of gyration and
        first-to-last distance that follow from them, and finite extensibility
        (no pair further apart than its contour length). The expected values are
        those of the discretized chain ``sample_conformations()`` draws from,
        which differs from the continuous worm-like chain by at most 0.02% (the
        report notes the exact figure). The model's analytical end-to-end
        distribution (O'Brien 2009) and the Benoit-Doty radius of gyration are reported as context.

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
        m, correlation = self._ensemble_discretization()
        discrete = mean_squared_distance_map(self.nres, lambda k: discrete_worm_like_chain_msd(k, self.b, m, correlation))
        continuous = self.get_mean_squared_distance_map()
        off_diagonal = continuous > 0
        deviation = float(np.max(np.abs(discrete[off_diagonal] / continuous[off_diagonal] - 1))) if np.any(off_diagonal) else 0.0

        # the O'Brien P(r) can always be evaluated for this chain (the constructor
        # has already checked it is at least one persistence length long), but not
        # necessarily for the N-1 residues between the end beads, so that context
        # comparison is skipped when it cannot be evaluated
        reference_re = float(self.get_root_mean_squared_end_to_end_distance())
        end_to_end = None
        if self.nres >= 2:
            try:
                end_to_end = WormLikeChain2('A' * (self.nres - 1), self.p_of_r_resolution, lp=self.lp,
                                        aa_size=self.b).get_end_to_end_distribution()
            except WLC2Exception:
                pass
        reference_rg = float(self.get_mean_radius_of_gyration())

        expectations = ModelExpectations(
            "Worm-like chain (O'Brien)", discrete,
            contour_length_per_residue=self.b,
            reference_rg=reference_rg,
            reference_rg_label='RMS Rg (Benoit-Doty, N residues)',
            reference_re=reference_re,
            reference_re_label='whole-chain RMS Re from its analytical P(r) (N residues)',
            end_to_end_distribution=end_to_end,
            end_to_end_label="the O'Brien 2009 P(r) for N-1 residues (an approximation to the worm-like chain)",
            notes=(f'The ensemble is a discretized worm-like chain ({m} straight sub-segments per residue); its mean-squared '
                   f'distances differ from the continuous worm-like chain by at most {100 * deviation:.3f}%, and the checks above '
                   f'use the discretized values.',))
        return compare_ensemble_to_model(conformations, expectations)

