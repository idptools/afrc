"""
afrc.py
An analytical version of the Flory Random Coil (FRC) for polypeptides, implemented using the rotational isomeric state approximation of Flory and Volkenstein and parameterized on the excluded volumed dihedral backbone maps.

"""
import numpy as np
from .polymer import PolymerObject
from .config import P_OF_R_RESOLUTION, AA_list
from .iofunctions import validate_keyword


class AFRCException(Exception):
    """
    Exception class specific for the Analytical FRC

    """
    pass

            
class AnalyticalFRC:
    """
    The AnalyticalFRC object is the only user-facing object that the AFRC package provides. All functionality
    is associated with function called from this object, and the object itself is instantiated with a single
    amino acid string. For all intents and purposes, one can think of as an *AnalyticalFRC* object as holding
    one protein sequence and providing an interface to ask specific types of polymer questions.

    .. code-block:: python

          from afrc import AnalyticalFRC
          MyProtein = AnalyticalFRC('KFGGPRDQGSRHDSEQDNSDNNTIFVQGLG')
    
    Note
    ----
    Distributions and parameters are only calculate as requested, such that initializing an Analytical FRC 
    object is a cheap operation. However, operations relating to intramolecular distances (``get_distance_map()``, 
    ``get_internal_scaling()`` etc. are more computationally expensive.


    Attributes
    ----------    
    seq : str
       The amino acid sequence of the protein being examined. 

        
    """



    # .....................................................................................
    #
    def __init__(self, seq, adaptable_P_res=False):
        """
        Constructor for an AFRC object which can be queried to obtain varies parameters and statistics.

        Parameters
        ----------
        seq : str
            Amino acid sequence for the protein of interest (case insensitive). If this is an invalid string 
            it will raise an AFRCException.

        adaptable_P_res : Bool (False)
           Sets the resolution used for generating probability distributions. By default this is assigned to
           a fixed value (0.05 A). However, if this is set to True a sequence-specific adaptable resolution
           is used and calculated as :math:`d_{max} / 500.00` (where :math:`d_{max}` reflects the contour 
           length of the polypeptide and is defined as :math:`3.7n`.
           
        
        """

        try:
            seq = seq.upper()
        except AttributeError as e:
            raise AFRCException('Input must be a string of amino acids')
            
        # check a valid string was passed and assign to object variable 
        self.__check_seq_is_valid(seq)            
        self.seq = seq 

        # set up what our P of R spacing will look like...
        if adaptable_P_res:
            dmax=3.7*len(seq)
            self.p_of_r_resolution = dmax/500.0
        else:
            self.p_of_r_resolution = P_OF_R_RESOLUTION


        self.full_seq_PO = PolymerObject(seq, self.p_of_r_resolution)
        self.matrix=False



    # .....................................................................................
    #
    def __check_seq_is_valid(self, seq):
        """
        Helper function that ensures the passed sequence contains ONLY valid amino acids (upper-case).

        No return value and totally stateless.

        """
        for i in seq: 
            if i not in AA_list:
                raise AFRCException('Passed amino acid sequence contains non-standard amino acids [%s]' %(i))



    # .....................................................................................
    #
    def __build_matrix(self):
        """

        Internal function that limits matrix construction until its actually needed! The matrix in question
        here is an [n x n] matrix of PolymerObjects for querying inter-residue distances.

        """

        if self.matrix is False:
            self.matrix = []
            for i in range(0, len(self.seq)):
                row = []
                for j in range(0, len(self.seq)):
                    row.append(0)
                self.matrix.append(row)
            
            for i in range(0, len(self.seq)):
                for j in range(i, len(self.seq)):
                    subseq = self.seq[i:j]                
                    self.matrix[i][j] = PolymerObject(subseq, self.p_of_r_resolution)
                    self.matrix[j][i] = self.matrix[i][j]
        else:
            pass



    # .....................................................................................
    #        
    def __len__(self):
        """
        Returns the length of the sequence
        """
        return len(self.seq)



    # .....................................................................................
    #
    def get_distance_map(self, calculation_mode='scaling law'):
        """
        Returns the complete inter-residue distance map, an [n x n] upper-right triangle
        matrix that can be used as a reference set for constructing scaling maps.

        Distances are in angstroms and are measured from the residue center of mass.
        
        Parameters
        ----------


        Returns
        -------
        np.ndarray
           An [n x n] square matrix (where n = length of the amino acid sequence) defining
           the inter-residue distances between every pair of residues. 
        
        """

        # check input mode information
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        # construct the internal matrix of polymers
        self.__build_matrix()

        # initialize the distance-distance matrix
        dm = np.zeros((len(self.seq),len(self.seq)))
        
        # for each inter-residue dustance (not only upper right
        # triangle is computed)
        for i in range(0,len(self.seq)):
            for j in range(i,len(self.seq)):
                dm[i,j] = self.matrix[i][j].get_mean_end_to_end_distance(calculation_mode)

        return dm



    # .....................................................................................
    #
    def get_internal_scaling(self, calculation_mode='scaling law'):
        """
        Returns the internal scaling profile - a [2 by n] matrix that reports on the average
        distance between all residues that are n positions apart ( where n  is | i - j | ). 

        Distances are in angstroms and are measured from the residue center of mass.

        A linear log-log fit of this data gives a gradient of 0.5 (:math:`\\nu^{app} = 0.5`).
        
        Parameters
        ----------

        Returns
        -------
        np.ndarray
           An [2 x n] matrix (where n = length of the amino acid sequence). The first column
           is the set of | i-j | distances, and the second defines the average inter-residue 
           distances between every pair of residues that are | i-j | residues apart in sequnce 
           space.                
        
        """

        # validate mode and construct the matrix if not yet built
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')
        self.__build_matrix()

        # set the empty dictionary and iterate through all non-redundant distances        
        rij={}


        # now cycle through every non-redundant pair
        for i in range(0,len(self.seq)):
            for j in range(i+1,len(self.seq)):
                
                # if empty initialize 
                if j-i not in rij:
                    rij[j-i] = []

                rij[j-i].append(self.matrix[i][j].get_mean_end_to_end_distance(calculation_mode))
                    
        # having established all possible distances we then
        # calculate the average 
        k = list(rij)
        k.sort()
        mean_vals= []
        for dis in k:
            mean_vals.append(np.mean(rij[dis]))

        return np.array((k,mean_vals)).transpose()
     
                   

    # .....................................................................................
    #            
    def get_rg_distribution(self):        
        """
        Defines the radius of gyration (:math:`R_g`) distribution using equation (3) from [Lhuillier1988]_. 

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        return self.full_seq_PO.get_radius_of_gyration_distribution()



    # .....................................................................................
    #
    def get_re_distribution(self):        
        """
        Defines the end-to-end distance (Re) distribution using the standard end-to-end model (as in [Rubinstein2003]_). 
        
        :math:`P(r)=4\pi r^2 \Biggl( \\frac {3} {2 \pi \\langle r^2 \\rangle} \Biggr)^2 e^{\\frac{3r^2}{2\\langle r^2 \\rangle}}`


        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and 
           the second array is the probability of that distance.

        """

        return self.full_seq_PO.get_end_to_end_distribution()



    # .....................................................................................
    #
    def get_re_distribution_WLC(self):        
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

        return self.full_seq_PO.get_end_to_end_distribution_WLC()



    # .....................................................................................
    #        
    def get_mean_rg(self, calculation_mode='scaling law'):
        """
        Returns the mean radius of gyration (:math:`R_g`) as calculated from the 
        :math:`R_g` distribution.

        Paramaters
        ...........


        calculation_mode : string (keyword)
             calculation_mode defines the mode in which the average is calculated, and can be 
             set to either 'scaling law' (default) or 'distribution'. If 'distribution' is used
             then the complete Rg distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Rg = R0 * N^{0.5} is used.        

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_radius_of_gyration(calculation_mode)


    def get_mean_re(self, calculation_mode='scaling law'):
        """
        Returns the mean end-to-end distance (:math:`R_e`). This value can be the absolute
        mean end-to-end distance or the root-mean-sequence end-to-end distance.

        Parameters
        ----------

        calculation_mode : string (keyword)
             calculation_mode defines the mode in which the average is calculated, and can be 
             set to either 'scaling law' (default) or 'distribution'. If 'distribution' is used
             then the complete Re distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Re = R0 * N^{0.5} is used.        


        Returns
        -------
        float
           Value equal to the average end-to-end distance (as defined by ``mode``).

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_end_to_end_distance(calculation_mode)

    # .....................................................................................
    #
    def get_mean_re_WLC(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`). As calculated from the Worm-like
        chain (WLC) model as defined by Zhou [Zhou2004]_. 
        

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        return self.full_seq_PO.get_mean_end_to_end_distance_WLC()
        

    # .....................................................................................
    #
    def get_interresidue_distance_distribution(self, R1, R2):
        """
        Returns the distribution between a pair of residues on the chain.

        
        Parameters
        ----------

        R1 : int
           The first residue of the pair being investigated.

        R2 : int
           The second residue of the pair being investigated.


        Returns
        -------

        np.ndarray 
           Probability distribution returned as a 2D numpy array in which column one is 
           the distances (in Angstroms) and colum two is the probablility.
           
        """                
        self.__build_matrix()
        return self.matrix[R1][R2].get_end_to_end_distribution()



    # .....................................................................................
    #
    def get_mean_interresidue_distance(self, R1, R2, calculation_mode='scaling law'):

        """
        Returns the mean distance between a pair of residues on the chain.
        
        Parameters
        ----------

        R1 : int
           The first residue of the pair being investigated.

        R2 : int
           The second residue of the pair being investigated.

        calculation_mode : string (keyword)
             calculation_mode defines the mode in which the average is calculated, and can be 
             set to either 'scaling law' (default) or 'distribution'. If 'distribution' is used
             then the complete Rg distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Rg = R0 * N^{0.5} is used.        


        Returns
        -------

        float
           Absolute mean distance (if mode='mean_distance') or root mean-squared 
           
           
        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.matrix[R1][R2].get_mean_end_to_end_distance(calculation_mode)


    # .....................................................................................
    #

    def get_mean_interresidue_radius_of_gyration(self, R1, R2, calculation_mode='scaling law'):
        """
        Returns the mean radius of gyration (:math:`R_g`) as calculated from the 
        :math:`R_g` distribution BETWEEN a pair of residues (i.e. the :math:`R_g` distribution
        for an internal local region of the chain).

        Parameters
        ----------

        calculation_mode : string (keyword)
             calculation_mode defines the mode in which the average is calculated, and can be 
             set to either 'scaling law' (default) or 'distribution'. If 'distribution' is used
             then the complete Rg distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Rg = R0 * N^{0.5} is used.        


        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.matrix[R1][R2].get_mean_radius_of_gyration(calculation_mode)
        

    # .....................................................................................
    #
    def sample_rg_distribution(self,n=1000):
        """
        Subsamples from the :math:`R_g` distirbution to generate an uncorrelated 'trajectory' 
        of points. Useful for creating a sized-match sample to compare with simulation
        data.

        Parameters
        ----------
        n : int
           Number of random values to sample (default = 1000)

        Returns
        -------
        np.ndarray
           Returns an n-length array with n independent values (floats)

        """

        return self.full_seq_PO.sample_radius_of_gyration_distribution(dist_size=n)



    # .....................................................................................
    #
    def sample_re_distribution(self,n=1000):
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

        return self.full_seq_PO.sample_end_to_end_distribution(dist_size=n)




