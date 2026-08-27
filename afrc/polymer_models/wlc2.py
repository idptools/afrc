import numpy as np
from afrc.config import P_OF_R_RESOLUTION

class WLC2Exception(Exception):
    pass

class WormLikeChain2:
    """
    This class generates an object that returns polymer statistics consistent with the Worm-like chain
    model as implemented by O'Brien. Provides mean Re, mean Rg, and Re distribution.


    [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). 
    How accurate are polymer models in the analysis of Forster resonance 
    energy transfer experiments on proteins? The Journal of Chemical Physics, 
    130(12), 124903.

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

        # set sequence info
        self.nres = len(seq)

        # first cast to floats
        self.lp = float(lp)
        self.b = float(aa_size)

        # also sanity check

        if self.lp <= 0:
            raise WLC2Exception('Error, lp cannot be less than or equal to 0')

        if self.b <= 0:
            raise WLC2Exception('Error, aa_size cannot be less than or equal to 0')

        Lc = self.b*self.nres

        # the chain must be at least one persistence length long. Note this compares
        # the CONTOUR length (N*aa_size, in Angstroms) with the persistence length -
        # this check previously compared the number of residues against lp, which mixes
        # a residue count with a length in Angstroms
        if Lc < self.lp:
            raise WLC2Exception('Passed sequence has a contour length (%.2f A) shorter than the persistence length (%.2f A)' % (Lc, self.lp))

        # next calculate params as defined by O'Brien et al
        
        self.alpha = (3*Lc) / (4*self.lp)
        
        self.C2 = 1/(2*self.lp)

        t1 = np.power(np.pi, 3/2)
        t2 = np.exp(-self.alpha)
        t3 = np.power(self.alpha,-3/2)

        # warning - if 
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
        O'Brien et al.

        This is a composition independent model for which the end-to-end distance depends
        solely on the number of amino acids. It is included here as an additional reference 
        model.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        References
        -----------
        [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009). 
        How accurate are polymer models in the analysis of Forster resonance 
        energy transfer experiments on proteins? The Journal of Chemical Physics, 
        130(12), 124903.


        """
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by O'brien et al.

        Note, the mean here is calculated by integrating over P(r) vs r.
        
        Returns
        -------
        float
           Value equal to the mean end-to-end distance distribution

        """
        [a,b] = self.get_end_to_end_distribution()

        return np.sum(a*b)


    # .....................................................................................
    #        
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by O'brien et al.

        Note mean here is calculated by taking the square root after integrating over P(r) vs r^2.
        
        Returns
        -------
        float
           Value equal to the root-mean-squared end-to-end distance

        """


        [a,b] = self.get_end_to_end_distribution()
        return np.sqrt(np.sum(b*np.power(a,2)))



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
        prefactor = np.min((np.max([4,Lp]),10))*0.5

        p_dist = np.arange(0, prefactor*(7*np.power(self.nres,0.5)), self.p_of_r_resolution)

        # precompute the prefactor
        PREFACT = np.pi*self.C1*4

        # compute P(r) across the whole grid at once. Beyond the contour length
        # (1 - (r/Lc)^2) is negative and the expression is undefined (nan) - the
        # chain cannot be longer than its contour length, so those points carry zero
        # probability; the errstate silences the associated warnings
        r2 = np.power(p_dist,2)
        RoL2 = np.power(p_dist/Lc,2)
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            LHS = (PREFACT*r2) / (Lc*np.power((1-RoL2),9/2))
            RHS = (-3*Lc) / (4*Lp*(1-RoL2))
            p_val_raw = LHS*np.exp(RHS)

        p_val_raw = np.nan_to_num(p_val_raw, nan=0.0, posinf=0.0, neginf=0.0)
        p_val_raw[p_dist >= Lc] = 0.0

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist


    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self):
        """
        Returns the mean radius of gyration (:math:`R_g`) as defined by
        O'Brien et al in [1]. NOTE it doesn't explicitly say it in the
        paper, but we're assuming this is actually Rg^{2} so this returns
        the square root of the Rg defined in table 1 (WLC row).

        With :math:`C_2 = 1/(2L_p)` the expression used is

            <Rg^2> = Lc/(6 C2) - 1/(4 C2^2) + 1/(4 C2^3 Lc)
                     - (1 - exp(-Lc/Lp)) / (8 C2^4 Lc^2)

        which is the standard Benoit-Doty worm-like chain result
        :math:`\\langle R_g^2 \\rangle = L_c L_p/3 - L_p^2 + 2L_p^3/L_c - 2L_p^4/L_c^2 (1 - e^{-L_c/L_p})`
        rewritten in terms of :math:`C_2`. In the rigid-rod limit
        (:math:`L_c \\ll L_p`) it correctly reduces to :math:`L_c^2/12`.

        [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009).
        How accurate are polymer models in the analysis of Forster resonance
        energy transfer experiments on proteins? The Journal of Chemical Physics,
        130(12), 124903.

        [2] Benoit, H., & Doty, P. (1953). Light scattering from non-Gaussian
        chains. The Journal of Physical Chemistry, 57(9), 958-963.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        Lc = self.nres*self.b
        C2 = self.C2
        Lp = self.lp

        # NOTE: the second term carries a MINUS sign (it is -Lp^2 in the
        # Benoit-Doty form); it was previously written as +1/(4*C2^2), which
        # over-estimated Rg (badly so for short chains, where it left a
        # spurious constant 2*Lp^2 offset instead of the correct Lc^2/12
        # rigid-rod limit)
        return np.sqrt(Lc/(6*C2) - 1/(4*np.power(C2,2)) +  1/(Lc*4*np.power(C2, 3)) - (1 - np.exp(-Lc/Lp))/(8*np.power(C2, 4)*np.power(Lc, 2)))
