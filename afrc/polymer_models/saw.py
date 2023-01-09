import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class SAWException:
    pass

class SAW:
    """
    This class generates an object that returns polymer statistics consistent with a 
    self-avoiding random walk (SAW). This model was developed by Jhulian 'J' Alston,
    and is based on the reference implementation by O'Brein et al [1]

    [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). 
    How accurate are polymer models in the analysis of Forster resonance 
    energy transfer experiments on proteins? The Journal of Chemical Physics, 
    130(12), 124903.


    """

    # .....................................................................................
    #        
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION):
        """
        Method to create SAW (self-avoiding walk) object. Seq should be a valid upper-case 
        amino acid sequence and p_of_r_resolution defines the resolution (in angstroms) to 
        be used for distributions.
        
        By default p_of_r_resolution is taken from the config.py file in the afrc package which 
        defines the resolution at 0.05 A.
        
        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate number of residues)

        p_of_r_resolution : float
            Bin width for bulding probability distributions. In Angstroms.


        """

        # constants derived by 
        self.a = 3.67853
        self.b = 1.23152

        self.theta = 0.3
        self.delta = 2.5

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
    def get_end_to_end_distribution(self, prefactor=5.5):
        """
        Defines the end-to-end distribution based on the SAW as defined by
        https://aip.scitation.org/doi/10.1063/1.3082151. 

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. It is included here as an additional reference 
        model.

        By default this uses a prefactor of 5.5 A (0.55 nanometers).

        Parameters
        ------------
        prefactor : float
            Prefactor is a number that tunes the SAW dimensions. 0.5 is in the right ballpark
            but this number should be tuned to match EV sims.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        self.__compute_end_to_end_distribution(prefactor)

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self, prefactor=5.5):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated 
        from the SAW model as defined 
        https://aip.scitation.org/doi/10.1063/1.3082151. 

        By default this uses a prefactor of 5.5 A (0.55 nanometers).
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        
        [a,b] = self.get_end_to_end_distribution(prefactor)
        
        return np.sqrt(np.sum(np.power(a,2)* b))

    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self, prefactor=5.5):
        """
        
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        gamma = 1.1615
        nu=0.589
        top = gamma*(gamma + 1)
        bottom = 2*(gamma + 2*nu)*(gamma + 2*nu + 1)

        Ree = self.get_mean_end_to_end_distance(prefactor=prefactor)

        return np.sqrt(Ree**2*(top/bottom))


    
    # .....................................................................................
    #        
    def __compute_end_to_end_distribution(self, prefactor):
        """
        Defines the end-to-end distribution based on the SAW as defined by
        https://aip.scitation.org/doi/10.1063/1.3082151
        . This is where we actually perform the polymer model calculation.

        
        """
    
        # define the chainlength-dependent prefactor
        Ree = prefactor*np.power(self.nres, 0.598)

        # use same pdist as was used for the parent AFRC model - this is the set of (r) values
        # on a P(r) vs. r plot
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # initialize an empty array
        p_val_raw = np.zeros(len(p_dist))
        
        # define SAW normalization factors as defined by 
        # https://aip.scitation.org/doi/10.1063/1.3082151
        a = self.a
        b = self.b

        theta = self.theta
        delta = self.delta

        # for each possible value of 'r'
        for i in range(0,len(p_dist)):

            r = p_dist[i] 
            
            #compute p(r)
            P_r_one = a/Ree
            P_r_two = np.power((r/Ree), theta+2)
            P_r_three = np.exp(-b*(np.power((r/Ree), delta)))
            
            p_val_raw[i] = (P_r_one)*(P_r_two)*(P_r_three)
            #print(p_val_raw)
            
        p_val_length = len(p_val_raw)

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist[:p_val_length]
