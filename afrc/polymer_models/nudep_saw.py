import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice
from scipy.special import gamma as GAMMA_FUNCTION

class NuDepSAWException(Exception):
    pass

class NuDepSAW:
    """
    This class generates an object that returns polymer statistics consistent with a
    nu-dependent self-avoiding random walk (SAW), as developed by Zheng et al. and
    written in the form used here by Soranno.

    This is the same universal scaling form used by the fixed-exponent
    :class:`~afrc.polymer_models.saw.SAW` model, but with the Flory scaling exponent
    (nu) left as a free parameter, so a single model spans the range from a collapsed
    globule (nu ~ 1/3), through the ideal/theta chain (nu = 0.5), to a fully solvated
    good-solvent coil (nu ~ 0.588).

    This is a composition-independent model: the sequence is used only to set the number
    of residues. It is included as an additional reference model.

    [1] Zheng, W., Zerze, G. H., Borgia, A., Mittal, J., Schuler, B., & Best, R. B.
    (2018). Inferring properties of disordered chains from FRET transfer efficiencies.
    The Journal of Chemical Physics, 148(12), 123329.

    [2] Soranno, A. (2020). Physical basis of the disorder-order transition. Archives of
    Biochemistry and Biophysics, 685, 108305.

    """

    # .....................................................................................
    #        
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION):
        """
        Method to create nu-dependent SAW (self-avoiding walk) object. Seq should be a valid 
        upper-case  amino acid sequence and p_of_r_resolution defines the resolution (in 
        angstroms) to  be used for distributions.
        
        By default p_of_r_resolution is taken from the config.py file in the afrc package which 
        defines the resolution at 0.05 A.
        
        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate number of residues)

        p_of_r_resolution : float
            Bin width for bulding probability distributions. In Angstroms.

        Returns
        -----------
        Generates an NuDepSAW object


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
        Internal function for computing the first prefactor term  (A1) in equation 9a/9b
        in Soranno (2020). Delta and g are defined as:

        delta = 1/(1-nu)
        g = (gamma-1)/nu

        A1 is the normalization constant that makes P(r) integrate to one, and is

            A1 = (delta/4pi) * Gamma[(5+g)/delta]^((3+g)/2) / Gamma[(3+g)/delta]^((5+g)/2)

        Note the gamma-function arguments are (5+g)/delta and (3+g)/delta, and the
        (3+g)/2 and (5+g)/2 terms are EXPONENTS rather than multipliers.

        Parameters
        -------------
        delta : float
            First parameters

        g : float
            second parameters

        Returns
        ---------------
        Float
            Returns the A1 prefactor


        References
        -------------
        Soranno, A. (2020). Physical basis of the disorder-order transition. Archives of
        Biochemistry and Biophysics, 685, 108305.
        """

        T1 = delta/(4*np.pi)
        T2_top = np.power(GAMMA_FUNCTION((5+g)/delta), (3+g)/2)
        T2_bottom = np.power(GAMMA_FUNCTION((3+g)/delta), (5+g)/2)

        return T1* (T2_top/T2_bottom)


    # .....................................................................................
    #            
    def __compute_A2(self, delta, g):
        """
        Internal function for computing the second prefactor term (A2) in equation 9a/9b
        in Soranno (2020). Delta and g are defined as:

        delta = 1/(1-nu)
        g = (gamma-1)/nu

        A2 sets the width of the distribution and is fixed by requiring that the
        root-mean-square end-to-end distance of P(r) is exactly Ree, giving

            A2 = ( Gamma[(5+g)/delta] / Gamma[(3+g)/delta] )^(delta/2)

        As for A1, the gamma-function arguments are (5+g)/delta and (3+g)/delta.

        Parameters
        -------------
        delta : float
            First parameters

        g : float
            second parameters

        Returns
        ---------------
        Float
            Returns the A2 prefactor

        References
        -------------
        Soranno, A. (2020). Physical basis of the disorder-order transition. Archives of
        Biochemistry and Biophysics, 685, 108305.
        """

        top    = GAMMA_FUNCTION((5+g)/delta)
        bottom = GAMMA_FUNCTION((3+g)/delta)

        return np.power(top/bottom, delta/2)


    # .....................................................................................
    #        
    def get_end_to_end_distribution(self, nu=0.5, prefactor=5.5):
        """
        Defines the end-to-end distribution based on the nu-dependent SAW model.

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. Both nu and the prefactor can be varied, and
        together they set the root-mean-square end-to-end distance to
        ``prefactor * N^nu``.

        Parameters
        ------------
        nu : float
            Flory scaling exponent. Must lie strictly between 0 and 1; physically
            meaningful values fall between 0.33 and 0.6.

        prefactor : float
            Prefactor is a number that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        # note this model does not memoize because nu and prefactor can change
        # so we don't 
        self.__compute_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return (self.__p_of_Re_R, self.__p_of_Re_P)

    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self, nu=0.5, prefactor=5.5):
        """
        Returns the mean end-to-end distance (:math:`R_e`) for the nu-dependent
        SAW model.

        The mean is computed by integrating over the :math:`P(r)` vs. :math:`r`
        distribution (i.e. :math:`\\sum r \\cdot P(r)`), consistent with the
        convention used by the other models in this package.

        Parameters
        ----------
        nu : float
            Flory scaling exponent. Should fall between 0.33 and 0.6.

        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the mean end-to-end distance.

        """

        [a, b] = self.get_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return np.sum(a * b)

    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self, nu=0.5, prefactor=5.5):
        """
        Returns the root-mean-square end-to-end distance
        (:math:`\\sqrt{\\langle R_e^2 \\rangle}`) for the nu-dependent SAW model.

        The value is computed by taking the square root after integrating over
        :math:`P(r)` vs. :math:`r^2`.

        Parameters
        ----------
        nu : float
            Flory scaling exponent. Should fall between 0.33 and 0.6.

        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the root-mean-square end-to-end distance.

        """

        [a, b] = self.get_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return np.sqrt(np.sum(np.power(a, 2) * b))

    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self, nu=0.5, prefactor=5.5):
        """
        Returns the mean radius of gyration (:math:`R_g`) for the nu-dependent
        SAW model.

        :math:`R_g` is obtained from the mean-squared end-to-end distance via the
        analytical ratio :math:`\\langle R_g^2 \\rangle / \\langle R_e^2 \\rangle`
        expressed in terms of the gamma exponent and the scaling exponent
        :math:`\\nu`.

        Parameters
        ----------
        nu : float
            Flory scaling exponent. Should fall between 0.33 and 0.6.

        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

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
        Subsamples from the end-to-end distance distribution to generate an uncorrelated 
        'trajectory' of points. Useful for creating a size-matched sample to compare with
        simulation data.

        Parameters
        ----------
        n : int
           Number of random values to sample (default = 1000)

        nu : float
            Flory scaling exponent. Should fall between 0.33 and 0.6

        prefactor : float
            Prefactor is a number that tunes the SAW dimensions. Default is 5.5 A.

        
        Returns
        -------
        np.ndarray
           Returns an n-length array with n independent values (floats)

        """

        # note this model does not memoize because nu and prefactor can change
        # so we don't 
        self.__compute_end_to_end_distribution(nu=nu, prefactor=prefactor)

                
        return choice(self.__p_of_Re_R, n, p=self.__p_of_Re_P)


    

    # .....................................................................................
    #        
    def __compute_end_to_end_distribution(self, nu=0.5, prefactor=5.5):
        """
        Defines the end-to-end distribution based on the nu-dependent SAW 
        as defined by Zheng et al.

        Parameters
        -----------------
        nu : float
            Flory scaling exponent

        Prefactor : float
            Prefactor 
        
        

        This is where we actually perform the polymer model calculation.


        The specific expression used here comes from the equation formulated
        by Soranno (Eq 9b in [2])


        References
        ---------------
        [1] 
        

        [2] Soranno, A. (2020). Physical basis of the disorder-order transition.
        Archives of Biochemistry and Biophysics, 685, 108305.


        """

        # nu must sit strictly between 0 and 1 - at nu = 1 delta blows up and at
        # nu = 0 g blows up, so guard rather than emit a ZeroDivisionError
        nu = float(nu)
        if nu <= 0 or nu >= 1:
            raise NuDepSAWException('Error, nu must be between 0 and 1 (physically meaningful values run from ~0.33 to ~0.6)')

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
