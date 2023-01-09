import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice
from scipy.special import gamma as GAMMA_FUNCTION

class NuDepSAWException:
    pass

class NuDepSAW:
    """
    This class generates an object that returns polymer statistics consistent with a 
    nu-dependent self-avoiding random walk (SAW). As developed by Zheng et al 2020


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

        # set gamma - orginally defined in
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
        T2_top = GAMMA_FUNCTION(5+(g/delta))*((3+g)/2)
        T2_bottom = GAMMA_FUNCTION(3+(g/delta))*((5+g)/2)

        return T1* (T2_top/T2_bottom)
    

    # .....................................................................................
    #            
    def __compute_A2(self, delta, g):
        """
        Internal function for computing the second prefactor term (A2) in equation 9a/9b
        in Soranno (2020). Delta and g are defined as:

        delta = 1/(1-nu)
        g = (gamma-1)/nu

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

        top    = GAMMA_FUNCTION(5+(g/delta))
        bottom = GAMMA_FUNCTION(3+(g/delta))

        return np.power(top/bottom, delta/2)    
            
    # .....................................................................................
    #        
    def get_end_to_end_distribution(self, nu=0.5, prefactor=5.5):
        """
        Defines the end-to-end distribution based 

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. Both nu and the prefactor can be varied 
        

        Parameters
        ------------
        nu : float
            Flory scaling exponent. Should fall between 0.33 and 0.6

        prefactor : float
            Prefactor is a number that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        # note this model does not memoize necause nu and prefactor can change
        # so we don't 
        self.__compute_end_to_end_distribution(nu=nu, prefactor=prefactor)

        return (self.__p_of_Re_R, self.__p_of_Re_P)

    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self, nu=0.5, prefactor=5.5):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the SAW model as defined 
        https://aip.scitation.org/doi/10.1063/1.3082151. 
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        
        [a,b] = self.get_end_to_end_distribution(nu=nu, prefactor=prefactor)
        
        return np.sqrt(np.sum(np.power(a,2)* b))


    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self, nu=0.5, prefactor=5.5):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the SAW model as defined 
        https://aip.scitation.org/doi/10.1063/1.3082151. 
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        top = self.gamma*(self.gamma + 1)
        bottom = 2*(self.gamma + 2*nu)*(self.gamma + 2*nu + 1)

        Ree = self.get_mean_end_to_end_distance(nu=nu, prefactor=prefactor)

        return np.sqrt(Ree**2*(top/bottom))
    

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
        
        gamma = self.gamma
        g = (gamma -1)/nu
        delta = 1/(1-nu)
        A1 = self.__compute_A1(delta, g)
        A2 = self.__compute_A2(delta, g)            
    
        # define the chainlength-dependent prefactor
        ##
        ## WARNING: For reasons I cannot explain, one needs to multiply the prefactor
        ## by pi to get the distance units right. I suspect this reflects an error
        ## somewhere, but multiplying by pi gives us quantitative agreement with both
        ## SAW (when prefactors are the same and nu=0.598, within error of the approximations
        ## made, AND with the AFRC when nu=0.5
        ##
        ## If anyone figures out where this mysterious pi factor has come from please let
        ## me know!
        ##
        ## ~ash Jan 2023
        Ree = prefactor*np.power(self.nres, nu)*np.pi

        # use same pdist as was used for the parent AFRC model - this is the set of (r) values
        # on a P(r) vs. r plot
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # initialize an empty array
        p_val_raw = np.zeros(len(p_dist))

        # first term in EQ 9b as written by Soranno 2020)
        T1 = (A1*4*np.pi)/Ree


        # for each possible value of 'r'
        for i in range(0,len(p_dist)):

            r = p_dist[i]

            # second term in EQ 9b (as written by Soranno 2020)
            T2 = np.power(r/Ree, 2+g)

            # third term in EQ 9b (as written by Soranno 2020)
            T3 = np.exp(-A2*np.power(r/Ree,delta))
            
            
            p_val_raw[i] = T1*T2*T3

            
        p_val_length = len(p_val_raw)

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist[:p_val_length]
