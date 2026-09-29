"""
polymer.py

The ``polymer`` module contains the ``PolymerObject`` class, which holds the
AFRC description of a single chain (or chain segment) of fixed sequence.

``AnalyticalFRC`` builds one ``PolymerObject`` for the full sequence and, when
inter-residue quantities are requested, one for every sub-sequence between a
pair of residues. Users do not normally need to create these directly.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from .config import AA_list, RIJ_RMS_R0, RIJ_R0, RG_X0, RG_R0, P_OF_R_RESOLUTION
from .exceptions import AFRCException
from numpy.random import choice

class PolymerObject:
    """
    Internal object that describes one chain (or chain segment) under the AFRC.

    On construction we compute the composition-weighted prefactors for the
    sequence. Distributions are only built when they are requested, and are then
    cached unless ``cache_distributions=False``.

    A zero-length sequence is allowed (the inter-residue matrix needs one for
    each residue paired with itself). All of its distributions put their weight
    at zero.

    Attributes
    ----------
    nres : int
        Number of residues in the sequence.

    zero_length : bool
        True if the sequence is empty.

    RMS_Re_scaling : float
        The root-mean-square end-to-end distance, :math:`R_0^{rms} N^{0.5}`,
        in Angstroms (0 for a zero-length chain).

    p_of_r_resolution : float
        Grid spacing (in Angstroms) used for the distributions.

    cache_distributions : bool
        Whether distributions are kept once computed.

    """

    # .....................................................................................
    #
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION, cache_distributions=True):
        """
        Create a PolymerObject for a sequence.

        Parameters
        ----------
        seq : str
            Valid upper-case amino acid sequence. May be empty. Note this is not
            validated here - ``AnalyticalFRC`` does that.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distributions. Default is
            0.05 A (from ``config.py``).

        cache_distributions : bool
            If True (default) each distribution is kept after it is first
            computed, so later calls are free. If False it is recomputed on every
            call and never stored. ``AnalyticalFRC`` turns caching off for the
            [n x n] inter-residue matrix: recomputing one pair's distribution
            takes tens of microseconds, whereas keeping all of them costs
            gigabytes of memory for a few hundred residues.

        """

        # set sequence info
        self.nres = len(seq)
        self.zero_length = False
        self.RMS_Re_scaling = 0
        self.p_of_r_resolution = p_of_r_resolution
        self.cache_distributions = cache_distributions

        # set distribution info to false - is calculated as needed
        self.__p_of_Re_R = False
        self.__p_of_Re_P = False

        self.__p_of_Rg_R = False
        self.__p_of_Rg_P = False


        ## *********************************
        ## Construct sequence-specific prefactors
        ##

        # inter-residue distance prefactors
        self.__R0_RMS = 0
        self.__R0 = 0

        # RG prefactor for Lhuillier equation 'X0'
        self.__X0 = 0

        # RG prefactor for the Rg scaling law (<Rg> = RG_R0 * N^{0.5})
        self.__RG_R0 = 0

        # if sequence is empty no need to compute anything, and set the
        # zero length flag to to true, and return (we're done!). This allows
        # the code to natively deal with zero-length strings rather than throwing
        # an exception.
        if len(seq) == 0:
            self.zero_length = True
            return

        # Compute the sequence-specific prefactors using the global lookup tables. This is
        # just calculating the compositionally-weighted average value for the R0_RMS, R0 and X0
        # prefactors. Recall that
        #
        #  Root mean squared Re = R0_RMS*N^{0.5}
        #  <Re> = R0*N^{0.5}
        #  <Rg> = RG_R0*N^{0.5}
        #
        # (X0 is the prefactor that enters the Lhuillier P(Rg) expression rather than
        # a scaling-law prefactor, see __compute_Rg_distribution)
        #
        # Note we count each amino acid once and skip any that are absent. This
        # matters because the inter-residue matrix builds ~n^2/2 of these objects -
        # counting every residue type four times over made building the matrix
        # O(n^3). Skipping absent residues only drops terms that add exactly 0.0,
        # and the loop order is unchanged, so the prefactors are bitwise identical
        for AA in AA_list:
            count = seq.count(AA)
            if count == 0:
                continue

            fraction = count/float(self.nres)
            self.__R0_RMS = self.__R0_RMS + fraction*RIJ_RMS_R0[AA]
            self.__R0 = self.__R0 + fraction*RIJ_R0[AA]
            self.__RG_R0 = self.__RG_R0 + fraction*RG_R0[AA]

            # note - we apply a +0.005 offset to each RG_X0 value
            self.__X0 = self.__X0 + fraction*(RG_X0[AA]+0.005)

        # and then compute the ensemble average RMS-Re and the absolute ensemble average Re
        # using the standard scaling law (R0 * N^{nu}) where nu=0.5 and R0 is calculated
        # based on composition
        self.RMS_Re_scaling = self.__R0_RMS * np.power(self.nres,0.5)

    # .....................................................................................
    #
    def get_end_to_end_distribution(self):
        """
        Return the Gaussian-chain end-to-end distance distribution.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1).

        """

        # return the cached distribution if we have one, otherwise compute it (the
        # compute function caches it if this object is caching)
        if self.__p_of_Re_R is not False:
            return (self.__p_of_Re_R, self.__p_of_Re_P)

        return self.__compute_end_to_end_distribution()



    # .....................................................................................
    #
    def get_radius_of_gyration_distribution(self):
        """
        Return the Lhuillier radius of gyration distribution.

        Returns
        -------
        tuple of np.ndarray
            ``(radii, probabilities)``, where radii are in Angstroms and the
            probabilities are a normalized probability mass function (they sum
            to 1).

        """

        # return the cached distribution if we have one, otherwise compute it (the
        # compute function caches it if this object is caching)
        if self.__p_of_Rg_R is not False:
            return (self.__p_of_Rg_R, self.__p_of_Rg_P)

        return self.__compute_Rg_distribution()



    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self, calculation_mode='scaling law'):
        """
        Return the mean end-to-end distance.

        Parameters
        ----------
        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle R_e \\rangle = R_0 N^{0.5}`, or ``'distribution'``,
            which takes the expectation over the end-to-end distribution.

        Returns
        -------
        float
            The mean end-to-end distance (in Angstroms).

        Raises
        ------
        AFRCException
            If an unrecognized ``calculation_mode`` is passed.

        """

        # if we're using the scaling law relationship
        if calculation_mode == 'scaling law':
            return self.__R0 * np.power(self.nres,0.5)

        # if we're calculating the expected value from the distribution
        elif calculation_mode == 'distribution':
            [a,b] = self.get_end_to_end_distribution()
            return np.sum(a*b)

        else:
            raise AFRCException(f"calculation_mode must be set to one of ['distribution', 'scaling law'] (was set to {calculation_mode})")




    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self, calculation_mode='scaling law'):
        """
        Return the mean radius of gyration.

        Parameters
        ----------
        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle R_g \\rangle = R_0^{g} N^{0.5}` with the
            composition-weighted ``RG_R0`` prefactor, or ``'distribution'``, which
            takes the expectation over the radius of gyration distribution.

        Returns
        -------
        float
            The mean radius of gyration (in Angstroms).

        Raises
        ------
        AFRCException
            If an unrecognized ``calculation_mode`` is passed.

        """

        # if we're using the scaling law relationship. Note we use the calibrated Rg
        # prefactor here rather than the ideal-chain relation Rg = <Re>/sqrt(6); the
        # latter relates root-mean-square (not mean) radii and so under-estimates <Rg>
        # by ~5%, which would leave the two calculation_modes disagreeing
        if calculation_mode == 'scaling law':
            return self.__RG_R0 * np.power(self.nres, 0.5)

        # if we're calculating the expected value from the distribution
        elif calculation_mode == 'distribution':
            (a,b) = self.get_radius_of_gyration_distribution()
            return np.sum(a*b)

        else:
            raise AFRCException(f"calculation_mode must be set to one of ['distribution', 'scaling law'] (was set to {calculation_mode})")



    # .....................................................................................
    #
    def sample_end_to_end_distribution(self, dist_size=1000):
        """
        Draw random end-to-end distances from the distribution.

        Parameters
        ----------
        dist_size : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``dist_size`` independent end-to-end distances (in Angstroms). For a
            zero-length chain every value is 0.

        """
        if self.zero_length:
            return np.repeat(0.0,dist_size)
        else:
            (r, p) = self.get_end_to_end_distribution()
            return choice(r, dist_size, p=p)

    # .....................................................................................
    #
    def sample_radius_of_gyration_distribution(self, dist_size=1000):
        """
        Draw random radii of gyration from the distribution.

        Parameters
        ----------
        dist_size : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``dist_size`` independent radii of gyration (in Angstroms). For a
            zero-length chain every value is 0.

        """
        if self.zero_length:
            return np.repeat(0.0,dist_size)
        else:
            (r, p) = self.get_radius_of_gyration_distribution()
            return choice(r, dist_size, p=p)



    # .....................................................................................
    #
    def __compute_end_to_end_distribution(self):
        """
        Build the Gaussian-chain end-to-end distribution (and cache it, if caching).

        .. math::

           P(r) = 4\\pi r^2 \\left( \\frac{3}{2\\pi \\langle R_e^2 \\rangle} \\right)^{3/2}
                  \\exp\\left( -\\frac{3 r^2}{2 \\langle R_e^2 \\rangle} \\right)

        with :math:`\\langle R_e^2 \\rangle = (R_0^{rms})^2 N`. The distribution is
        evaluated from 0 to :math:`21\\sqrt{N}` Angstroms (over three times the
        root-mean-square size) and normalized to sum to 1.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``.

        """

        # a zero-length chain has no extent, so all its weight sits at r = 0. Handling
        # this here (rather than letting the arithmetic below divide by zero) is what
        # lets zero-length PolymerObjects - which the [n x n] inter-residue matrix
        # necessarily contains along its diagonal - behave sanely
        if self.zero_length:
            return self.__store_end_to_end_distribution(np.array([0.0]), np.array([1.0]))

        # set distance range we're going to calculate P of R over = max is 3* peak - a somewhat
        # arbitrarily big, safe value. Note we employ a heuristic to ensure for short chains this
        # remains sufficient
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # compute the ensemble average square end-to-end distance. Note that
        # self.__R0_RMS * np.power(self.nres,0.5) gives SQRT(<Re^2>), so by squaring
        # this we get the correct parameter (i.e. mean-squared end-to-end distance)
        self.mean_squared_re = np.power(self.__R0_RMS * np.power(self.nres,0.5),2)

        # evaluate the Gaussian chain P(r) across the whole grid at once (this used to
        # go through np.vectorize, which made distribution-mode distance maps ~100x
        # slower than they need to be)
        A = np.power((3/(2*np.pi*self.mean_squared_re)),3.0/2.0)
        p_val_raw = 4*np.pi*A*np.power(p_dist,2)*np.exp(-(3*p_dist*p_dist)/(2*self.mean_squared_re))

        return self.__store_end_to_end_distribution(p_dist, p_val_raw/np.sum(p_val_raw))


    # .....................................................................................
    #
    def __store_end_to_end_distribution(self, distances, probabilities):
        """
        Cache an end-to-end distribution if this object is caching, and return it.

        Parameters
        ----------
        distances : np.ndarray
            Distance grid (in Angstroms).

        probabilities : np.ndarray
            Normalized probabilities on that grid.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, unchanged.

        """
        if self.cache_distributions:
            self.__p_of_Re_R = distances
            self.__p_of_Re_P = probabilities

        return (distances, probabilities)


    # .....................................................................................
    #
    def __store_Rg_distribution(self, radii, probabilities):
        """
        Cache a radius of gyration distribution if this object is caching, and return it.

        Parameters
        ----------
        radii : np.ndarray
            Radius of gyration grid (in Angstroms).

        probabilities : np.ndarray
            Normalized probabilities on that grid.

        Returns
        -------
        tuple of np.ndarray
            ``(radii, probabilities)``, unchanged.

        """
        if self.cache_distributions:
            self.__p_of_Rg_R = radii
            self.__p_of_Rg_P = probabilities

        return (radii, probabilities)


    # .....................................................................................
    #
    def __compute_Rg_distribution(self):
        """
        Build the radius of gyration distribution (and cache it, if caching).

        This uses equation 3 of Lhuillier (1988) with :math:`\\nu = 0.5` and
        :math:`d = 3`, evaluated at :math:`\\rho = X_0 R_g` where :math:`X_0` is
        the composition-weighted prefactor. The distribution is evaluated from 0
        to :math:`6\\sqrt{N}` Angstroms and normalized to sum to 1.

        Returns
        -------
        tuple of np.ndarray
            ``(radii, probabilities)``.

        References
        ----------
        [1] Lhuillier, D. (1988). A simple model for polymeric fractals in a good
        solvent and an improved version of the Flory approximation. Journal de
        Physique, 49(5), 705-710.

        """

        # as for the end-to-end distribution, a zero-length chain puts all its
        # weight at zero
        if self.zero_length:
            return self.__store_Rg_distribution(np.array([0.0]), np.array([1.0]))

        # note 0.5 reflects the scaling exponent, 3 is the dimensionality
        alpha = 1/(0.5*3 - 1)
        delta = 1/(1-0.5)

        # setup r values and the empty np vector we're gonna populate
        # p_val_rg_r = np.arange(0,self.nres, self.p_of_r_resolution) old way...
        p_val_rg_r = np.arange(0,3*(2*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # N raised to the power of nu (0.5)
        N_nu = np.power(self.nres,0.5)

        # evaluate the Lhuillier expression across the whole grid at once. The
        # errstate means we don't complain about the divide-by-zero at r_mod = 0,
        # where (N_nu/r_mod) is infinite, exp(-inf) is zero and P(0) is correctly 0
        r_mod = p_val_rg_r * self.__X0
        with np.errstate(divide='ignore'):
            f = np.exp(-np.power(N_nu/r_mod, alpha*3) - np.power((r_mod/N_nu), delta))
            p_val_raw = np.power(self.nres, -0.5*3)*f*(r_mod/N_nu)

        return self.__store_Rg_distribution(p_val_rg_r, p_val_raw/np.sum(p_val_raw))



    # .....................................................................................
    #
    def get_mean_inverse_end_to_end_distance(self):
        """
        Return the mean inverse end-to-end distance, :math:`\\langle 1/R_e \\rangle`.

        For the Gaussian end-to-end distribution this has the closed form

        .. math::

           \\langle 1/R_e \\rangle = \\sqrt{\\frac{6}{\\pi \\langle R_e^2 \\rangle}}

        which follows from integrating :math:`P(r)/r`. This is the quantity the
        Kirkwood-Riseman hydrodynamic radius needs. Note it is *not* the same as
        :math:`1/\\langle R_e \\rangle` - for a Gaussian chain the two differ by a
        factor of exactly :math:`4/\\pi`.

        Returns
        -------
        float
            The mean inverse end-to-end distance (in inverse Angstroms).

        Raises
        ------
        AFRCException
            If the polymer has zero length, for which the quantity is undefined.

        """

        if self.zero_length:
            raise AFRCException('The mean inverse distance is undefined for a zero-length polymer')

        # RMS_Re_scaling is sqrt(<Re^2>) from the composition-weighted scaling law
        return np.sqrt(6.0/(np.pi*self.RMS_Re_scaling*self.RMS_Re_scaling))


    # .....................................................................................
    #
    def compute_apparent_rms_bond_length(self):
        """
        Return the apparent root-mean-square bond length of the chain.

        For an ideal chain :math:`\\langle R_e^2 \\rangle = (N-1) b^2`, so this
        returns :math:`b = \\sqrt{\\langle R_e^2 \\rangle / (N-1)}`.

        Returns
        -------
        float
            The apparent root-mean-square bond length (in Angstroms).

        Raises
        ------
        AFRCException
            If the polymer has fewer than two residues (and so no bonds).

        """

        # with fewer than two residues there are zero bonds, so (N-1) is zero or
        # negative and the expression below is undefined - fail loudly rather than
        # returning an inf/nan
        if self.nres < 2:
            raise AFRCException('Cannot compute an apparent bond length for a polymer with fewer than two residues')

        re_v = np.sqrt((self.RMS_Re_scaling * self.RMS_Re_scaling) / (self.nres - 1))
        return re_v
