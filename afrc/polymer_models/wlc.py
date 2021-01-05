import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class WormLikeChain:
    """
    This class generates an object that returns polymer statistics consistent with the Worm-like chain
    model as implemented by Houx (2004)

    """

    # .....................................................................................
    #        
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION):
        """
        Method to create Polymer Object. Seq should be a valid upper-case amino acid sequence and p_of_r_resolution
        defines the resolution (in angstroms) to be used for distributions.

        By default p_of_r_resolution is taken from the config.py file in the afrc package which defines the resolution
        at 0.05 A.

        """

        # note that input validation is done in the AnalyticalFRC object constructor

        # set sequence info
        self.nres = len(seq)
        self.zero_length = False
        self.p_of_r_resolution = p_of_r_resolution

        # set distribution info to false - these are calculated if/when needed
        self.__p_of_Re_R = False
        self.__p_of_Re_P = False

        if len(seq) == 0:
            self.zero_length = True


    # .....................................................................................
    #        
    def get_end_to_end_distribution(self):

        """
        Defines the end-to-end distribution based on the Worm-like chain (WLC) as defined by
        Zhou [Zhou2004]_. 

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. It is included here as an additional reference 
        model.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        # if we have not yet computed the WLC end-to-end distance distribution do it now
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by Zhou [Zhou2004]_. 
        
        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        [a,b] = self.get_end_to_end_distribution()
        return np.sum(a*b)


    ##########################################################################################
    ##
    def __compute_end_to_end_distribution(self):
        """
        Defines the end-to-end distribution based on the Worm-like chain (WLC) as defined by
        Zhou. This is where we actually perform the polymer model calculation.

        """

        # define persistence length and contour length
        Lp = 3.0
        Lc = self.nres*3.8

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


