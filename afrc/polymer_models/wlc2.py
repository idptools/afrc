import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class WLC2Exception:
    pass

class WormLikeChain2:
    """
    This class generates an object that returns polymer statistics consistent with the Worm-like chain
    model as implemented by Houx (2004)

    """

    # .....................................................................................
    #        
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION, lp=3.0, aa_size=3.8):
        """
        Method to create Polymer Object. Seq should be a valid upper-case amino acid sequence and p_of_r_resolution
        defines the resolution (in angstroms) to be used for distributions.

        By default p_of_r_resolution is taken from the config.py file in the afrc package which defines the resolution
        at 0.05 A.

        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate number of residues)

        p_of_r_resolution : float
            Bin width for bulding probability distributions. In Angstroms.

        lp : float
            Persistence length. We use a default of 3 but 4 is also used a lot in the literature.

        aa_size : float
            Size of one amino acid (called 'b' in the literature). 3.8 is the generally acceptable 
            value used.

        """

        raise WLC2Exception("THIS IS NOT WORKING RIGHT BUT I DON'T KNOW WHY ARGHHHHHHHHHHH")

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
        
        Lc = self.b*self.nres
        # next calculate params as defined by O'Brien et al
        self.alpha = (3*Lc) / (4*self.lp)
        self.C2 = 1/(2*self.lp)

        t1 = np.power(np.pi, 3/2)
        t2 = np.exp(-self.alpha)
        t3 = np.power(self.alpha,-3/2)
        t4 = 3*np.power(self.alpha,-1)
        t5 = (15/4)*np.power(self.alpha,-2)

        self.C1 = np.power(t1*t2*t3*(1 + t4 + t5), -1)

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
        # by having the code like this it means we only ever compute the distribution 
        # once (wich, if the calculation is expensive is good)
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


    # .....................................................................................
    #        
    def __compute_end_to_end_distribution(self):
        """
        Defines the end-to-end distribution based on the Worm-like chain (WLC) as defined by
        O'Brien. This is where we actually perform the polymer model calculation.

        """

        # define persistence length and contour length
        Lp = self.lp
        Lc = self.nres*self.b
        

        # use same pdist as was used for the parent AFRC model 
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # initialize an empty array
        p_val_raw = np.zeros(len(p_dist))

        # precompute the prefactor
        PREFACT = np.pi*self.C1*4
        

        # for each possible value of 'r'
        for i in range(0,len(p_dist)):
            r = p_dist[i]            
            r2 = np.power(r,2)
            RoL2 = np.power(r/Lc,2)
            

            LHS = PREFACT*r2 /(Lc*np.power(1-RoL2,9/2))
            RHS = (-3*Lc) / (4*Lp*(1-RoL2))

            p_val_raw[i] = LHS*np.exp(RHS)
            
            
        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


    # .....................................................................................
    #        
    def get_mean_rg(self):
        """
        Returns the mean radius of gyration (:math:`R_g`) as defined by 
        O'brien et al in [1]


        [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). 
        How accurate are polymer models in the analysis of Forster resonance 
        energy transfer experiments on proteins? The Journal of Chemical Physics, 130(12), 124903.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        Lc = self.nres*self.b
                                                                                        

                                                                                              

        return Lc/(6*self.C2) + 1/(4*np.power(self.C2,2)) +  1/(Lc*4*np.power(self.C2,3)) - (1-np.exp(-Lc/self.lp))/(8*np.power(self.C2,4)*np.power(Lc,2))
