import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class FJCException(Exception):
    pass

class FreelyJointedChain:
    """
    This class generates an object that returns polymer statistics consistent with a
    freely jointed chain (FJC) of ``N`` rigid segments each of length ``b``.

    Unlike the (Gaussian) Analytical Flory Random Coil, the end-to-end distance
    distribution used here is the non-Gaussian Kuhn-Grün distribution [1], which
    is exact in the long-chain limit and - crucially - respects the finite
    extensibility of the chain (the end-to-end distance can never exceed the
    contour length :math:`L = Nb`). At small fractional extensions it reduces
    smoothly to the Gaussian result, so for typical IDP-length chains the bulk of
    the distribution is close to the AFRC, with deviations appearing in the tail.

    This is a composition-independent model: the sequence is used only to set the
    number of segments. It is included as an additional reference model.

    [1] Kuhn, W., & Grün, F. (1942). Beziehungen zwischen elastischen Konstanten
    und Dehnungsdoppelbrechung hochelastischer Stoffe. Kolloid-Zeitschrift,
    101(3), 248-271.

    [2] Cohen, A. (1991). A Padé approximant to the inverse Langevin function.
    Rheologica Acta, 30(3), 270-273.

    """

    # .....................................................................................
    #
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION, b=3.8):
        """
        Method to create a FreelyJointedChain object. Seq should be a valid upper-case
        amino acid sequence and p_of_r_resolution defines the resolution (in angstroms)
        to be used for distributions.

        By default p_of_r_resolution is taken from the config.py file in the afrc package
        which defines the resolution at 0.05 A.

        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate the number of segments).

        p_of_r_resolution : float
            Bin width for building probability distributions. In Angstroms.

        b : float
            Segment (Kuhn) length, in Angstroms. This is the FJC analogue of the
            ``aa_size`` parameter used by the worm-like chain models. The default of
            3.8 A corresponds to the Cα-Cα distance.

        """

        # cast and sanity check
        self.b = float(b)
        if self.b <= 0:
            raise FJCException('Error, b (segment length) cannot be less than or equal to 0')

        # set sequence info - the number of segments
        self.nres = len(seq)

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
        Defines the end-to-end distribution based on the freely jointed chain (FJC)
        using the non-Gaussian Kuhn-Grün distribution.

        This is a composition independent model for which the end-to-end distance
        depends solely on the number of amino acids. It is included here as an
        additional reference model.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and
           the second array is the probability of that distance.

        """
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`) for the freely jointed
        chain.

        The mean is computed by integrating over the :math:`P(r)` vs. :math:`r`
        distribution (i.e. :math:`\\sum r \\cdot P(r)`), consistent with the
        convention used by the other models in this package.

        Returns
        -------
        float
           Value equal to the mean end-to-end distance.

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sum(a * b)


    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Returns the root-mean-square end-to-end distance
        (:math:`\\sqrt{\\langle R_e^2 \\rangle}`) for the freely jointed chain.

        The value is computed by taking the square root after integrating over
        :math:`P(r)` vs. :math:`r^2`. In the long-chain (Gaussian) limit this
        approaches the ideal-chain result :math:`b\\sqrt{N}`.

        Returns
        -------
        float
           Value equal to the root-mean-square end-to-end distance.

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sqrt(np.sum(b * np.power(a, 2)))


    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self):
        """
        Returns the mean radius of gyration (:math:`R_g`) for the freely jointed
        chain.

        For an ideal chain the radius of gyration is related to the root-mean-square
        end-to-end distance by :math:`R_g = \\sqrt{\\langle R_e^2 \\rangle / 6}`, and
        this relationship is used here.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        return self.get_root_mean_squared_end_to_end_distance() / np.sqrt(6)


    # .....................................................................................
    #
    def sample_end_to_end_distribution(self, n=1000):
        """
        Subsamples from the end-to-end distance distribution to generate an uncorrelated
        'trajectory' of points. Useful for creating a sized-match sample to compare with
        simulation data.

        Parameters
        ----------
        n : int
           Number of random values to sample (default = 1000)

        Returns
        -------
        np.ndarray
           Returns an n-length array with n independent values (floats)

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
        Defines the end-to-end distribution based on the freely jointed chain (FJC).
        This is where we actually perform the polymer model calculation.

        The radial probability density is

            P(r) ∝ 4π r^2 exp( -N [ x β + ln(β / sinh β) ] )

        where x = r / (Nb) is the fractional extension and β = L^{-1}(x) is the
        inverse Langevin function, evaluated here using the Cohen Padé approximant
        β ≈ x(3 - x^2)/(1 - x^2). The distribution is only defined for r < Nb (the
        contour length), beyond which the probability is zero.

        """

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        # contour length
        L = self.nres * self.b

        # use the same style of r-grid as the other models, but never exceed the
        # contour length (the FJC distribution is undefined for r >= L)
        r_upper = min(3*(7*np.power(self.nres, 0.5)), L)
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
