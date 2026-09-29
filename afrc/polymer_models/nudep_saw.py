"""
nudep_saw.py

Self-avoiding walk model with a tunable Flory scaling exponent (nu).

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice
from scipy.special import gamma as GAMMA_FUNCTION

class NuDepSAWException(Exception):
    """Exception raised by the nu-dependent self-avoiding walk model."""
    pass

class NuDepSAW:
    """
    Self-avoiding walk with the Flory scaling exponent :math:`\\nu` as a free
    parameter, as developed by Zheng et al. [1] and written by Soranno [2].

    This is the same universal scaling form used by the fixed-exponent
    :class:`~afrc.polymer_models.saw.SAW`, but with :math:`\\nu` left free, so a
    single model spans a collapsed globule (:math:`\\nu \\approx 1/3`), the ideal
    chain (:math:`\\nu = 0.5`) and a good-solvent coil (:math:`\\nu \\approx 0.588`).
    Both ``nu`` and ``prefactor`` are passed to each method rather than set on the
    object.

    This is a composition-independent reference model: the sequence is only used
    to set the number of residues.

    References
    ----------
    [1] Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B.
    (2018). Inferring properties of disordered chains from FRET transfer
    efficiencies. The Journal of Chemical Physics, 148(12), 123329.

    [2] Soranno, A. (2020). Physical basis of the disorder-order transition.
    Archives of Biochemistry and Biophysics, 685, 108305.

    """

    # .....................................................................................
    #
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION):
        """
        Create a NuDepSAW object.

        Parameters
        ----------
        seq : str
            Amino acid sequence. Only its length is used, and it is not
            validated.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distribution. Default is
            0.05 A.

        """

        # set gamma - originally defined in
        # Le Guillou, J. C., & Zinn-Justin, J. (1977). Critical Exponents for the n-Vector
        # Model in Three Dimensions from Field Theory. Physical Review Letters, 39(2), 95–98.
        # for the case of n=0 (polymer), and raised in this context in the Soranno form
        # of the Zheng et al nu-dependent polymer model (see eq 9b in Soranno, A. (2020).
        # Physical basis of the disorder-order transition. Archives of Biochemistry and
        # Biophysics, 685, 108305.
        self.gamma = 1.1615

        # set sequence info
        self.nres = len(seq)

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
    def __compute_A1(self, delta, g):
        """
        Compute the normalization prefactor :math:`A_1` (Soranno 2020, Eq. 9b).

        .. math::

           A_1 = \\frac{\\delta}{4\\pi}
                 \\frac{\\Gamma[(5+g)/\\delta]^{(3+g)/2}}{\\Gamma[(3+g)/\\delta]^{(5+g)/2}}

        Note the gamma-function arguments are :math:`(5+g)/\\delta` and
        :math:`(3+g)/\\delta`, and :math:`(3+g)/2` and :math:`(5+g)/2` are
        exponents rather than multipliers.

        Parameters
        ----------
        delta : float
            Large-r exponent, :math:`\\delta = 1/(1 - \\nu)`.

        g : float
            Small-r exponent, :math:`g = (\\gamma - 1)/\\nu`.

        Returns
        -------
        float
            The :math:`A_1` prefactor, which normalizes :math:`P(r)` to 1.

        """

        T1 = delta/(4*np.pi)
        T2_top = np.power(GAMMA_FUNCTION((5+g)/delta), (3+g)/2)
        T2_bottom = np.power(GAMMA_FUNCTION((3+g)/delta), (5+g)/2)

        return T1* (T2_top/T2_bottom)


    # .....................................................................................
    #
    def __compute_A2(self, delta, g):
        """
        Compute the width prefactor :math:`A_2` (Soranno 2020, Eq. 9b).

        .. math::

           A_2 = \\left( \\frac{\\Gamma[(5+g)/\\delta]}{\\Gamma[(3+g)/\\delta]} \\right)^{\\delta/2}

        This value is fixed by requiring that the root-mean-square end-to-end
        distance of :math:`P(r)` is exactly :math:`R_{ee}`.

        Parameters
        ----------
        delta : float
            Large-r exponent, :math:`\\delta = 1/(1 - \\nu)`.

        g : float
            Small-r exponent, :math:`g = (\\gamma - 1)/\\nu`.

        Returns
        -------
        float
            The :math:`A_2` prefactor.

        """

        top    = GAMMA_FUNCTION((5+g)/delta)
        bottom = GAMMA_FUNCTION((3+g)/delta)

        return np.power(top/bottom, delta/2)


    # .....................................................................................
    #
    def get_end_to_end_distribution(self, nu=0.5, prefactor=5.5):
        """
        Return the end-to-end distance distribution.

        Because ``nu`` and ``prefactor`` can change between calls, the
        distribution is recomputed every time rather than cached.

        Parameters
        ----------
        nu : float
            Flory scaling exponent. Must lie strictly between 0 and 1;
            physically meaningful values run from ~0.33 to ~0.6. Default is 0.5.

        prefactor : float
            Size scale in Angstroms, such that the root-mean-square end-to-end
            distance is ``prefactor * N**nu``. Default is 5.5 A. Must be > 0.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1).

        Raises
        ------
        NuDepSAWException
            If ``nu`` is not strictly between 0 and 1, or ``prefactor`` is not
            positive.

        """

        # note this model does not memoize because nu and prefactor can change
        # so we don't
        self.__compute_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return (self.__p_of_Re_R, self.__p_of_Re_P)

    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self, nu=0.5, prefactor=5.5):
        """
        Return the mean end-to-end distance, :math:`\\langle R_e \\rangle`.

        This is the expectation over the end-to-end distribution,
        :math:`\\sum r P(r)`.

        Parameters
        ----------
        nu : float
            Flory scaling exponent, strictly between 0 and 1. Default is 0.5.

        prefactor : float
            Size scale in Angstroms. Default is 5.5 A.

        Returns
        -------
        float
            The mean end-to-end distance (in Angstroms).

        Raises
        ------
        NuDepSAWException
            If ``nu`` or ``prefactor`` is out of range.

        """

        [a, b] = self.get_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return np.sum(a * b)

    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self, nu=0.5, prefactor=5.5):
        """
        Return the root-mean-square end-to-end distance,
        :math:`\\sqrt{\\langle R_e^2 \\rangle}`.

        This is the square root of :math:`\\sum r^2 P(r)` over the end-to-end
        distribution, and equals ``prefactor * N**nu``.

        Parameters
        ----------
        nu : float
            Flory scaling exponent, strictly between 0 and 1. Default is 0.5.

        prefactor : float
            Size scale in Angstroms. Default is 5.5 A.

        Returns
        -------
        float
            The root-mean-square end-to-end distance (in Angstroms).

        Raises
        ------
        NuDepSAWException
            If ``nu`` or ``prefactor`` is out of range.

        """

        [a, b] = self.get_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return np.sqrt(np.sum(np.power(a, 2) * b))

    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self, nu=0.5, prefactor=5.5):
        """
        Return the root-mean-square radius of gyration,
        :math:`\\sqrt{\\langle R_g^2 \\rangle}`.

        This is obtained from the root-mean-square end-to-end distance via the
        universal ratio

        .. math::

           \\frac{\\langle R_g^2 \\rangle}{\\langle R_e^2 \\rangle} =
              \\frac{\\gamma(\\gamma + 1)}{2(\\gamma + 2\\nu)(\\gamma + 2\\nu + 1)}

        evaluated at the chosen :math:`\\nu`. Note that despite the method name
        this is the *root-mean-square* radius of gyration, not
        :math:`\\langle R_g \\rangle`. The name is kept for consistency with the
        other models.

        Parameters
        ----------
        nu : float
            Flory scaling exponent, strictly between 0 and 1. Default is 0.5.

        prefactor : float
            Size scale in Angstroms. Default is 5.5 A.

        Returns
        -------
        float
            The root-mean-square radius of gyration (in Angstroms).

        Raises
        ------
        NuDepSAWException
            If ``nu`` or ``prefactor`` is out of range.

        """

        top = self.gamma*(self.gamma + 1)
        bottom = 2*(self.gamma + 2*nu)*(self.gamma + 2*nu + 1)

        # the ratio above relates mean-squared radii, so we use the
        # root-mean-square end-to-end distance (sqrt(<Re^2>)) here
        Ree = self.get_root_mean_squared_end_to_end_distance(nu=nu, prefactor=prefactor)

        return np.sqrt(Ree**2*(top/bottom))

    # .....................................................................................
    #
    def sample_end_to_end_distribution(self, n=1000, nu=0.5, prefactor=5.5):
        """
        Draw random end-to-end distances from the distribution.

        Useful for building a size-matched sample to compare against simulation
        data.

        Parameters
        ----------
        n : int
            Number of values to draw. Default is 1000.

        nu : float
            Flory scaling exponent, strictly between 0 and 1. Default is 0.5.

        prefactor : float
            Size scale in Angstroms. Default is 5.5 A.

        Returns
        -------
        np.ndarray
            ``n`` independent end-to-end distances (in Angstroms). For a
            zero-length chain every value is 0.

        Raises
        ------
        NuDepSAWException
            If ``nu`` or ``prefactor`` is out of range.

        """

        # note this model does not memoize because nu and prefactor can change
        # so we don't
        self.__compute_end_to_end_distribution(nu=nu, prefactor=prefactor)


        return choice(self.__p_of_Re_R, n, p=self.__p_of_Re_P)




    # .....................................................................................
    #
    def __compute_end_to_end_distribution(self, nu=0.5, prefactor=5.5):
        """
        Build the end-to-end distribution for a given nu and prefactor.

        This is Eq. 9b of Soranno (2020):

        .. math::

           P(r) = \\frac{4\\pi A_1}{R_{ee}} \\left( \\frac{r}{R_{ee}} \\right)^{2+g}
                  \\exp\\left[ -A_2 \\left( \\frac{r}{R_{ee}} \\right)^{\\delta} \\right],
           \\qquad R_{ee} = \\texttt{prefactor}\\,N^{\\nu}

        with :math:`g = (\\gamma - 1)/\\nu` and :math:`\\delta = 1/(1 - \\nu)`.
        The grid runs from 0 to the larger of :math:`21\\sqrt{N}` and
        :math:`4R_{ee}`.

        Parameters
        ----------
        nu : float
            Flory scaling exponent, strictly between 0 and 1.

        prefactor : float
            Size scale in Angstroms. Must be > 0.

        Raises
        ------
        NuDepSAWException
            If ``nu`` or ``prefactor`` is out of range.

        """

        # nu must sit strictly between 0 and 1 - at nu = 1 delta blows up and at
        # nu = 0 g blows up, so guard rather than emit a ZeroDivisionError
        nu = float(nu)
        if nu <= 0 or nu >= 1:
            raise NuDepSAWException('Error, nu must be between 0 and 1 (physically meaningful values run from ~0.33 to ~0.6)')

        # the prefactor is a length scale, so it must be positive - a zero or negative
        # value previously produced an all-NaN distribution without complaint
        prefactor = float(prefactor)
        if prefactor <= 0:
            raise NuDepSAWException('Error, prefactor must be greater than 0 (it sets the size scale in Angstroms)')

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        gamma = self.gamma
        g = (gamma -1)/nu
        delta = 1/(1-nu)
        A1 = self.__compute_A1(delta, g)
        A2 = self.__compute_A2(delta, g)

        # define the chainlength-dependent prefactor. With A1/A2 computed correctly
        # this is exactly the root-mean-square end-to-end distance of the resulting
        # distribution (A2 is defined by that requirement), so no additional fudge
        # factor is needed here
        Ree = prefactor*np.power(self.nres, nu)

        # r values on a P(r) vs. r plot. We use the same grid as the parent AFRC model,
        # but ensure it always extends to at least 4*Ree - for long chains and/or large
        # nu, Ree grows faster than the sqrt(N) AFRC grid, and truncating the grid below
        # the tail of the distribution biases the mean and RMS values low
        upper = max(3*(7*np.power(self.nres, 0.5)), 4*Ree)
        p_dist = np.arange(0, upper, self.p_of_r_resolution)

        # first term in EQ 9b as written by Soranno 2020)
        T1 = (A1*4*np.pi)/Ree

        # second term in EQ 9b (as written by Soranno 2020)
        T2 = np.power(p_dist/Ree, 2+g)

        # third term in EQ 9b (as written by Soranno 2020)
        T3 = np.exp(-A2*np.power(p_dist/Ree, delta))

        p_val_raw = T1*T2*T3

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist
