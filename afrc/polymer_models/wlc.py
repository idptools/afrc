import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class WLCException:
    pass

class WormLikeChain:
    """
    This class generates an object that returns polymer statistics consistent with the Worm-like chain
    model as implemented by Zhou (2004).

    This model should be basically identical to the O'Brien model (WormLikeChain2) but show better
    numerical stability at large contour lengths. Unlike the O'Brien model this model does not
    provide an estimation of the mean Rg.

    Zhou, H.-X. (2004). Polymer models of protein stability, folding, and interactions. 
    Biochemistry, 43(8), 2141–2154.

    """

    # .....................................................................................
    #        
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION, lp=3.0, aa_size=3.8):
        """
        Method to create a WormLikeChain Object. Seq should be a valid upper-case amino acid 
        sequence and p_of_r_resolution defines the resolution (in angstroms) to be used for 
        distributions.
        
        By default p_of_r_resolution is taken from the config.py file in the afrc package 
        which defines the resolution at 0.05 A.
        

        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate number of residues)

        p_of_r_resolution : float
            Bin width for bulding probability distributions. In Angstroms.

        lp : float
            Persistence length. We use a default of 3 but 4 is also used a lot in the 
            literature.

        aa_size : float
            Size of one amino acid (called 'b' in the literature). 3.8 is the generally
            acceptable value used.
            

        """

        # note that input validation is done in the AnalyticalFRC object constructor

        # set sequence info
        self.nres = len(seq)

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
        Defines the end-to-end distribution based on the Worm-like chain (WLC).

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. It is included here as an additional reference 
        model.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by Zhou [Zhou2004]_. 

        Note mean here is calculated by integrating over P(r) vs r.
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        [a,b] = self.get_end_to_end_distribution()

        return np.sum(a*b)


    # .....................................................................................
    #        
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by Zhou [Zhou2004]_. 

        Note mean here is calculated by taking the square root after integrating over P(r) vs r^2.
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """


        [a,b] = self.get_end_to_end_distribution()
        return np.sqrt(np.sum(b*np.power(a,2)))

    # .....................................................................................
    #        
    def __compute_end_to_end_distribution(self):
        """
        Defines the end-to-end distribution based on the Worm-like chain (WLC) as defined by
        Zhou. This is where we actually perform the polymer model calculation.

        """

        # define persistence length and contour length
        Lp = self.lp
        Lc = self.nres*self.b

        # use same pdist as was used for the parent AFRC model 
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # initialize an empty array
        p_val_raw = np.zeros(len(p_dist))

        # precompute the prefactor
        prefactor_A = 4*np.pi*np.power(3.0/(4*np.pi*Lp*Lc),1.5)
        

        # define a function that depends on 'r'
        def zeta(r):
            return (1 - ((5*Lp/4*Lc) - 
                         ((2*np.power(r,2))/(np.power(Lc,2))) +
                         ((33*np.power(r,4))/(80*Lp*np.power(Lc,3))) +
                         ((79*np.power(Lp,2))/(160*np.power(Lc,2))) +
                         ((329*Lp*np.power(r,2))/(120*np.power(Lc,3))) -
                         ((6799*np.power(r,4))/(1600*np.power(Lc,4))) +
                         ((3441*np.power(r,6))/(2800*Lp*np.power(Lc,5))) -
                         ((1089*np.power(r,8))/(12800*np.power(Lp,2)*np.power(Lc,6)))))
            
            
        # for each possible value of 'r'
        for i in range(0,len(p_dist)):
            r = p_dist[i]            

            # compute P(r) at (r) based on the equations 5a/b in Zhou et al 2004
            p_val_raw[i] = prefactor_A*np.power(r,2)*np.exp(-3.0*(np.power(r,2))/(4*Lp*Lc))*zeta(r)
            
            
        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


