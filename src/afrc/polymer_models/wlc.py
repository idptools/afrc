"""
wlc.py

Worm-like chain (WLC) model using the closed-form approximation of Zhou (2004).

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from numpy.typing import NDArray
from afrc.ensemble import (discrete_worm_like_chain_msd, mean_squared_distance_map, sample_worm_like_chain, save_conformations,
                           validate_n_conformations, worm_like_chain_discretization, worm_like_chain_msd)
from afrc.ensemble_report import EnsembleReport, ModelExpectations, compare_ensemble_to_model
from afrc.config import P_OF_R_RESOLUTION

class WLCException(Exception):
    """Exception raised by the Zhou worm-like chain model."""
    pass

class WormLikeChain:
    """
    Worm-like chain model, as implemented by Zhou (2004).

    This is a composition-independent reference model: the sequence is only used
    to set the number of residues, and hence the contour length
    :math:`L_c = N b`. It should agree closely with the O'Brien model
    (``WormLikeChain2``), but unlike that model it does not provide a radius of
    gyration.

    The underlying expression is a series expansion in :math:`L_p/L_c` and
    :math:`r/L_c`, so it is only accurate when the contour length comfortably
    exceeds the persistence length (for the default parameters, chains of more
    than ~10-20 residues). No probability is assigned beyond the contour length,
    and a chain too short for the expansion to have any valid region raises a
    ``WLCException`` when the distribution is requested.

    References
    ----------
    [1] Zhou, H.-X. (2004). Polymer models of protein stability, folding, and
    interactions. Biochemistry, 43(8), 2141-2154.

    """

    # .....................................................................................
    #
    def __init__(self, seq: str, p_of_r_resolution: float = P_OF_R_RESOLUTION, lp: float = 3.0, aa_size: float = 3.8) -> None:
        """
        Create a WormLikeChain object.

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
        WLCException
            If ``lp`` or ``aa_size`` is not positive.

        """

        # note that input validation is done in the AnalyticalFRC object constructor

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
            raise WLCException('Error, lp cannot be less than or equal to 0')

        if self.b <= 0:
            raise WLCException('Error, aa_size cannot be less than or equal to 0')

        # p_of_r_resolution defines the P(r) resolution in angstroms - i.e. basically
        # the spacing between r values in a P(r) vs. r plot
        self.p_of_r_resolution = p_of_r_resolution

        # set distribution info to false - these are calculated if/when needed. M
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

        Raises
        ------
        WLCException
            If the contour length is too short relative to the persistence
            length for the Zhou expansion to be valid anywhere.

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
        Build and cache the end-to-end distribution (equations 5a/5b of Zhou 2004).

        .. math::

           P(r) = 4\\pi A r^2 \\exp\\left( -\\frac{3 r^2}{4 L_p L_c} \\right) \\zeta(r),
           \\qquad A = \\left( \\frac{3}{4\\pi L_p L_c} \\right)^{3/2}

        where :math:`\\zeta(r)` is Zhou's polynomial correction series. Points
        beyond the contour length and negative values of the series are set to
        zero, and the result is normalized to sum to 1.

        Raises
        ------
        WLCException
            If no part of the distribution survives (the contour length is
            comparable to or shorter than the persistence length).

        """

        # a zero-length chain has zero contour length (which would divide by zero
        # below), so put all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        # define persistence length and contour length
        Lp = self.lp
        Lc = self.nres*self.b

        # use the same r-grid as the parent AFRC model, but make sure it always extends
        # to at least four times the ideal-chain size scale sqrt(<r^2>) = sqrt(2*Lp*Lc).
        # The AFRC grid (21*sqrt(N)) is comfortably wide for the default lp = 3 A, but
        # for stiffer chains (lp of ~6 A and above) it cut into the tail of P(r) and
        # biased the mean and RMS values low
        upper = max(3*(7*np.power(self.nres, 0.5)), 4*np.sqrt(2*Lp*Lc))
        p_dist = np.arange(0, upper, self.p_of_r_resolution)

        # precompute the prefactor
        prefactor_A = 4*np.pi*np.power(3.0/(4*np.pi*Lp*Lc),1.5)

        # the polynomial correction series (equation 5b in Zhou 2004), evaluated across
        # the whole grid at once. Note the overall '1 - (...)' - the signs inside the
        # bracket are as written by Zhou, and this form reproduces the exact WLC
        # <r^2> = 2*Lp*Lc - 2*Lp^2*(1 - exp(-Lc/Lp)) to machine precision
        r = p_dist
        zeta = (1 - ((5*Lp/(4*Lc)) -
                     ((2*np.power(r,2))/(np.power(Lc,2))) +
                     ((33*np.power(r,4))/(80*Lp*np.power(Lc,3))) +
                     ((79*np.power(Lp,2))/(160*np.power(Lc,2))) +
                     ((329*Lp*np.power(r,2))/(120*np.power(Lc,3))) -
                     ((6799*np.power(r,4))/(1600*np.power(Lc,4))) +
                     ((3441*np.power(r,6))/(2800*Lp*np.power(Lc,5))) -
                     ((1089*np.power(r,8))/(12800*np.power(Lp,2)*np.power(Lc,6)))))

        # compute P(r) across the grid based on equations 5a/b in Zhou et al 2004
        p_val_raw = prefactor_A*np.power(r,2)*np.exp(-3.0*(np.power(r,2))/(4*Lp*Lc))*zeta

        # the end-to-end distance of a chain can never exceed its contour length. The
        # Zhou series is an expansion in Lp/Lc and r/Lc and is simply not meaningful
        # beyond r = Lc, where (for short and/or stiff chains) it can return substantial
        # positive values - for a 2-residue chain at lp = 3 A around 30% of the weight
        # sat beyond the contour length. Zero that region explicitly
        p_val_raw[r > Lc] = 0.0

        # the series can also produce spurious negative values in the far tail below Lc;
        # clamp these to zero so the result is a valid probability distribution
        p_val_raw[p_val_raw < 0] = 0.0

        # if nothing survives the two steps above the expansion has no valid region at
        # all - this happens when the contour length is comparable to or shorter than
        # the persistence length (a couple of residues at a large lp), which is outside
        # the regime the Zhou expansion is derived for. Fail loudly rather than
        # returning a NaN-filled distribution
        total = np.sum(p_val_raw)
        if not total > 0:
            raise WLCException('The Zhou (2004) worm-like chain expansion is not valid for a chain whose contour length (%.2f A) is comparable to or shorter than the persistence length (%.2f A)' % (Lc, Lp))

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/total
        self.__p_of_Re_R = p_dist


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
        distribution (Zhou 2004) is an approximation to it.

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
        distribution (Zhou 2004) is reported as context.

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

        # the analytical forms are only valid for long enough chains, so these
        # context comparisons are skipped when they cannot be evaluated
        reference_re = None
        end_to_end = None
        try:
            reference_re = float(self.get_root_mean_squared_end_to_end_distance())
        except WLCException:
            pass
        if self.nres >= 2:
            try:
                end_to_end = WormLikeChain('A' * (self.nres - 1), self.p_of_r_resolution, lp=self.lp,
                                        aa_size=self.b).get_end_to_end_distribution()
            except WLCException:
                pass
        reference_rg = None

        expectations = ModelExpectations(
            'Worm-like chain (Zhou)', discrete,
            contour_length_per_residue=self.b,
            reference_rg=reference_rg,
            reference_rg_label='RMS Rg (Benoit-Doty, N residues)',
            reference_re=reference_re,
            reference_re_label='whole-chain RMS Re from its analytical P(r) (N residues)',
            end_to_end_distribution=end_to_end,
            end_to_end_label='the Zhou 2004 P(r) for N-1 residues (an approximation to the worm-like chain)',
            notes=(f'The ensemble is a discretized worm-like chain ({m} straight sub-segments per residue); its mean-squared '
                   f'distances differ from the continuous worm-like chain by at most {100 * deviation:.3f}%, and the checks above '
                   f'use the discretized values.',))
        return compare_ensemble_to_model(conformations, expectations)

