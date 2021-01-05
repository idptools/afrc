"""polymer

The ``polymer`` module contains the ``PolymerObject`` class which contains information on a 
specific polymer of a fixed length. 



"""
import numpy as np
from .config import AA_list, RIJ_RMS_R0, RIJ_R0, RG_X0, RG_R0, P_OF_R_RESOLUTION
from numpy.random import choice

class PolymerObject:
    """
    Internal object 

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

        # set sequence info
        self.nres = len(seq)
        self.zero_length = False
        self.RMS_Re_scaling = 0
        self.p_of_r_resolution = p_of_r_resolution

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
        #  Root mean squared Re = N*R0_RMS^{0.5}
        #  <Re> = N*R0^{0.5}
        #  <Rg> = N*X0^{0.5}
        for AA in AA_list:
            self.__R0_RMS = self.__R0_RMS + (seq.count(AA)/float(self.nres))*RIJ_RMS_R0[AA]
            self.__R0 = self.__R0 + (seq.count(AA)/float(self.nres))*RIJ_R0[AA]

            # note - we apply a +0.005 offset to each RG_X0 value
            self.__X0 = self.__X0 + (seq.count(AA)/float(self.nres))*(RG_X0[AA]+0.005)

        # and then compute the ensemble average RMS-Re and the absolute ensemble average Re
        # using the standard scaling law (R0 * N^{nu}) where nu=0.5 and R0 is calculated
        # based on composition
        self.RMS_Re_scaling = self.__R0_RMS * np.power(self.nres,0.5)

    # .....................................................................................
    #        
    def get_end_to_end_distribution(self):
        """
        Function that returns the end-to-end distribution based on the analytical FRC
        model.

        Returns
        -------
        2D Numpy array in which the first column is the distance (in angstroms) and the second
        column is the probablity.
        
        
        """

        # if we have not yet done it, computed the end-to-end distance distribution
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)



    # .....................................................................................
    #        
    def get_radius_of_gyration_distribution(self):
        """
        Function that returns the radius of gyration distribution based on the analytical FRC
        model.
       

        Returns
        -------
        2D Numpy array in which the first column is the distance (in angstroms) and the second
        column is the probablity.

        """

        # if we have not yet computed the Rg distance distribution do so now
        if self.__p_of_Rg_R is False:
            self.__compute_Rg_distribution()

        return (self.__p_of_Rg_R, self.__p_of_Rg_P)



    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self, calculation_mode='scaling law'):
        """
        Function that returns the mean end-to-end distance using the analytical 
        FRC mode. Average can be computed using either the scalin law derivation
        or using the mean of the distribution.

        param
        

       

        Returns
        -------
        2D Numpy array in which the first column is the distance (in angstroms) and the second
        column is the probablity.


        """

        # if we're using the scaling law relationshup
        if calculation_mode == 'scaling law':
            return self.__R0 * np.power(self.nres,0.5)

        # if we're calculating the expected value from the distribution
        elif calculation_mode == 'distribution':
            [a,b] = self.get_end_to_end_distribution()
            return np.sum(a*b)




    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self, calculation_mode='scaling law'):
        """

        """

        # if we're using the scaling law relationshup
        if calculation_mode == 'scaling law':
            re = self.get_mean_end_to_end_distance('scaling law')
            return re/np.sqrt(6)
            
            
        # if we're calculating the expected value from the distribution
        elif calculation_mode == 'distribution':
            (a,b) = self.get_radius_of_gyration_distribution()
            return np.sum(a*b)

      

    # .....................................................................................
    #        
    def sample_end_to_end_distribution(self, dist_size=1000):
        """
        Function to randomly sample from the end-to-end distance distribution

        """
        if self.zero_length:
            return np.repeat(0.0,dist_size)
        else:
            if self.__p_of_Re_R is False:
                self.__compute_end_to_end_distribution()
                
            return choice(self.__p_of_Re_R, dist_size, p=self.__p_of_Re_P)
            
    ##########################################################################################
    ##
    def sample_radius_of_gyration_distribution(self, dist_size=1000):
        """
        Function to randomly sample from the Rg distance distribution

        """
        if self.zero_length:
            return np.repeat(0.0,dist_size)
        else:
            if self.__p_of_Rg_R is False:
                self.__compute_Rg_distribution()
            return choice(self.__p_of_Rg_R, dist_size, p=self.__p_of_Rg_P)



    ##########################################################################################
    ##
    def __compute_end_to_end_distribution(self):

        # set distance range we're going to calculate P of R over = max is 3* peak - a somewhat
        # arbitrarily big, safe value. Note we employ a heuristic to ensure for short chains this
        # remains sufficient
        p_dist = np.arange(0,3*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # initialize an empty vector which will be populated with probability 
        # values
        p_val_raw = np.zeros(len(p_dist))
        p_val_raw_alt = np.zeros(len(p_dist))
        
        # compute the ensemble average square end-to-end distance. Note that
        # self.__R0_RMS * np.power(self.nres,0.5) gives SQRT(<Re^2>), so by squaring
        # this we get the correct parameter (i.e. mean-squared end-to-end distance)
        
        self.mean_squared_re = np.power(self.__R0_RMS * np.power(self.nres,0.5),2)

        # precompute some stuff
        four_pi = 4*np.pi
        A = np.power((3/(2*np.pi*self.mean_squared_re)),3.0/2.0)

        def _prreturn(r):
            return four_pi*A*np.power(r,2)*np.exp(-(3*r*r)/(2*self.mean_squared_re))

        _prreturn_vectorized = np.vectorize(_prreturn)
        p_val_raw = _prreturn_vectorized(p_dist)
        
        # leave for now - old way - works but is a bit slower...
        """
        for i in xrange(0,len(p_dist)):
            #p_val[i] = 4*np.pi*np.power(p_dist[i],2)*np.power((3/(2*np.pi*mean_squared_re)),3.0/2)*np.exp(-(3*p_dist[i]*p_dist[i])/(2*mean_squared))
                p_val_raw_alt[i] = four_pi*A*np.power(p_dist[i],2)*np.exp(-(3*p_dist[i]*p_dist[i])/(2*self.mean_squared_re))
                if p_val_raw_alt[i] != p_val_raw[i]:
                    raise AFRCException('Diff in distributions -bug')
        """
                
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist
        
        
    ##########################################################################################
    ##
    def __compute_Rg_distribution(self):
        """
        Analytically computes the rg distribution


        Equation 3 from "A simple model for polymeric fractals in a good solvent and an improved version of the Flory approximation" by 
        Daniel Lhuillier, J. Phys. France 49 (1988) 705-710.

        """

        # note 0.5 reflects the scaling exponent, 3 is the dimensionality
        alpha = 1/((0.5*3 - 1));
        delta = 1/(1-0.5);
        
        # setup r values and the empty np vector we're gonna populate    
        # p_val_rg_r = np.arange(0,self.nres, self.p_of_r_resolution) old way...
        p_val_rg_r = np.arange(0,3*(2*np.power(self.nres,0.5)), self.p_of_r_resolution)
        p_val_raw = np.zeros(len(p_val_rg_r))

        # N raised to the power of nu (0.5)
        N_nu = np.power(self.nres,0.5)
        
        # define 
        def _prreturn(r_mod):            
            f = np.exp(   -   np.power(N_nu/r_mod, alpha*3) - np.power((r_mod / N_nu),delta))
            return (np.power(self.nres, -0.5*3)*f*(r_mod/N_nu))
            
        r_mod_vector = p_val_rg_r * self.__X0
        _prreturn_vectorized = np.vectorize(_prreturn)
        
        with np.errstate(divide='ignore'):
            p_val_raw_alt = _prreturn_vectorized(r_mod_vector)
            p_val_raw = p_val_raw_alt

        """
        # this with clause means we don't complain on div by zero when r_mod gets set to zero (i.e. at zero distance)
        with np.errstate(divide='ignore'):


            for i in xrange(0,len(p_val_rg_r)):

                r_mod = p_val_rg_r[i] * self.__X0 # this is where sequence specificty is applied
                f = np.exp(   -   np.power(N_nu/r_mod, alpha*3) - np.power((r_mod / N_nu),delta))

                p_val_raw[i] = np.power(self.nres, -0.5*3)*f*(r_mod/N_nu)                                
                #if p_val_raw[i] != p_val_raw_alt[i]:
                #    raise Exception('PROBLEMS')
        """


        self.__p_of_Rg_R = p_val_rg_r
        self.__p_of_Rg_P = p_val_raw/sum(p_val_raw)



    def compute_apparent_rms_bond_length(self):
        """

        """
        re_v =  np.sqrt((self.RMS_Re*self.RMS_Re)/(self.nres-1))
        return re_v

        
