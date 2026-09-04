"""
afrc.py
An analytical version of the Flory Random Coil (FRC) for polypeptides, implemented using the rotational isomeric state approximation of Flory and Volkenstein and parameterized on the excluded volume dihedral backbone maps.

"""
import numpy as np
from .polymer import PolymerObject
from .config import P_OF_R_RESOLUTION, AA_list
from .exceptions import AFRCException
from .iofunctions import validate_keyword
from .polymer_models import wlc

# AFRCException is defined in exceptions.py but re-exported here so that the
# long-standing `from afrc.afrc import AFRCException` import path keeps working
__all__ = ['AFRCException', 'AnalyticalFRC']


class AnalyticalFRC:
    """
    The AnalyticalFRC object is the main user-facing object that the AFRC 
    package provides. All functionality is associated with function called 
    from this object, and the object itself is instantiated with a single    
    amino acid string. For all intents and purposes, one can think of as 
    an *AnalyticalFRC* object as holding one protein sequence and providing
    an interface to ask specific types of polymer questions.

    .. code-block:: python

          from afrc import AnalyticalFRC
          MyProtein = AnalyticalFRC('KFGGPRDQGSRHDSEQDNSDNNTIFVQGLG')

    Note
    ----
    Distributions and parameters are only calculated as requested, such that
    initializing an AnalyticalFRC object is a cheap operation. However,
    operations relating to intramolecular distances (``get_distance_map()``,
    ``get_internal_scaling()`` etc.) are more computationally expensive.


        
    """


    # .....................................................................................
    #
    def __init__(self, seq, adaptable_P_res=False):
        """
        Constructor for an AFRC object which can be queried to obtain varies 
        parameters and statistics.

        Parameters
        ----------
        seq : str
            Amino acid sequence for the protein of interest (case insensitive).
            If this is an invalid (or empty) string it will raise an AFRCException.

        adaptable_P_res : Bool (False)
           Sets the resolution used for generating probability distributions. 
           By default this is assigned to a fixed value (0.05 A). However, if 
           this flag set to True a sequence-specific adaptable resolution           
           is used and calculated as :math:`d_{max} / 500.00` (where :math:`d_{max}` 
           reflects the contour length of the polypeptide and is defined as 
           :math:`3.7n`.
                              
        """

        try:
            seq = seq.upper()
        except AttributeError:
            raise AFRCException('Input must be a string of amino acids')

        # an empty sequence has no polymer to describe. Rather than let it through
        # (where most quantities silently come out as 0 and the Nygaard Rh as NaN)
        # reject it here. Note that zero-length PolymerObjects are still used
        # internally for the diagonal of the inter-residue matrix - that is
        # deliberate and unaffected by this check
        if len(seq) == 0:
            raise AFRCException('Input sequence must contain at least one amino acid')

        # check a valid string was passed and assign to object variable
        self.__check_seq_is_valid(seq)
        self.seq = seq 

        # set up what our P of R spacing will look like...
        if adaptable_P_res:
            dmax = 3.7*len(seq)
            self.p_of_r_resolution = dmax/500.0
        else:
            self.p_of_r_resolution = P_OF_R_RESOLUTION

        self.full_seq_PO = PolymerObject(seq, self.p_of_r_resolution)
        self.matrix=False

        # finally we define other polymer models which are attached as their own
        # class objects
        self.worm_like_chain = wlc.WormLikeChain(seq, self.p_of_r_resolution)




    # .....................................................................................
    #
    def __check_seq_is_valid(self, seq):
        """
        Helper function that ensures the passed sequence contains ONLY 
        valid amino acids (upper-case).

        No return value and totally stateless.

        Parameters
        ------------
        seq : string
            Amino acid sequence

        """
        for i in seq: 
            if i not in AA_list:
                raise AFRCException(f'Passed amino acid sequence contains non-standard amino acids [{i}]')



    # .....................................................................................
    #
    def __build_matrix(self):
        """

        Internal function that limits matrix construction until its actually needed! 
        The matrix in question here is an [n x n] matrix of PolymerObjects for querying 
        inter-residue distances. This is computationally a tad expensive to build, so this 
        function employs a memoization approach whereby IF the matrix is needed
        it is built, but only if

        """

        # if the matrix is not yet built
        if self.matrix is False:
            self.matrix = []

            ## the first set of for-loops initialize a matrix of lists
            # for each residue in the sequence
            for i in range(0, len(self.seq)):
                row = []

                # for each second residue in the sequence
                for j in range(0, len(self.seq)):
                    row.append(0)
                    
                self.matrix.append(row)

            ## the second set of for-loops defines the inter-residue
            ## distance for each unique pair of residues
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


    def __validate_residue_index(self, R):
        """
        Internal function that validates a passed residue index actually makes sense

        Parameters
        -----------
        R : int (or a type castable to integer)
            An integer to be used for residue selection

        Returns
        ----------
        int
            If this works returns the same integer (cast to an integer and everything!), else
            it will raises an exception 

        Raises
        ----------
        AFRCException

        """

        # note we catch TypeError as well as ValueError here - int('abc') raises the
        # former but int(None) raises the latter, and both should surface as an
        # AFRCException
        try:
            R = int(R)
        except (TypeError, ValueError):
            raise AFRCException('Could not convert residue [%s] to an integer' %(R))

        if R < 0:
            raise AFRCException('Residues %i cannot be under 0...'%(R))

        if R >= len(self):
            raise AFRCException('Residues %i cannot be over the chain length (%s)...'%(R, len(self)-1))

        return R
            


    # .....................................................................................
    #
    def get_distance_map(self, calculation_mode='scaling law', symmetric_map=False):
        """
        Returns the complete inter-residue distance map, an [n x n] upper-right triangle
        matrix that can be used as a reference set for constructing scaling maps.

        Distances are in angstroms and are measured from the residue center of mass.
        
        Parameters
        ----------
        calculation_mode : string (default = 'scaling law')
            A selector which must be equal to one of a specific set of options:

            'distribution' - means the P(r) distribution is used to calculate average distances

            'scaling law'  - means the derived scaling relationships are used to calculate the 
                             average distance

            If one of these is not provided then  an AFRCException is raised.

        symmetric_map : bool (default = False)
            If True, a full [n x n] matrix is returned, if False only the upper right triangle
            is returned.

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
        
        # for each inter-residue distance (only the upper-right
        # triangle is computed)
        for i in range(0,len(self.seq)):
            for j in range(i,len(self.seq)):
                dm[i,j] = self.matrix[i][j].get_mean_end_to_end_distance(calculation_mode)
                if symmetric_map:
                    dm[j,i] = dm[i,j]

        return dm



    # .....................................................................................
    #
    def get_internal_scaling(self, calculation_mode='scaling law'):
        """
        Returns the internal scaling profile - an [n-1 by 2] matrix that reports on the average
        distance between all residues that are k positions apart (where k is :math:`|i - j|`). 

        Distances are in angstroms and are measured from the residue center of mass.

        A linear log-log fit of this data gives a gradient of 0.5 (:math:`\\nu^{app} = 0.5`).

        Parameters
        ----------
        calculation_mode : string (keyword)
            calculation_mode defines the mode in which each inter-residue average is
            calculated, and can be set to either 'scaling law' (default) or
            'distribution'. If 'distribution' is used then the complete Re distribution
            is used to calculate the expected value. If the 'scaling law' is used then
            the standard Re = R0 * N^{0.5} is used.

        Returns
        -------
        np.ndarray
           An [n-1 x 2] matrix (where n = length of the amino acid sequence), with one row
           per sequence separation | i-j | from 1 to n-1. The first column is the set of
           | i-j | separations, and the second defines the average inter-residue distance
           between every pair of residues that are | i-j | residues apart in sequence space.
        
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
    def get_radius_of_gyration_distribution(self):        
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
    def get_end_to_end_distribution(self):
        r"""
        Defines the end-to-end distance (Re) distribution using the standard end-to-end model (as in [Rubinstein2003]_). 
        
        :math:`P(r) = 4\pi r^2 \left( \frac{3}{2\pi \langle r^2 \rangle} \right)^{3/2} e^{-\frac{3 r^2}{2 \langle r^2 \rangle}}`

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and
           the second array is the probability of that distance.

        """

        return self.full_seq_PO.get_end_to_end_distribution()



    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self, calculation_mode='distribution'):
        """
        Returns the mean radius of gyration (:math:`R_g`).

        Parameters
        ----------
        calculation_mode : str
             calculation_mode defines the mode in which the average is calculated, and can be
             set to either 'distribution' (default) or 'scaling law'. If 'distribution' is used
             then the complete Rg distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Rg = RG_R0 * N^{0.5} is used, where RG_R0
             is the composition-weighted radius of gyration prefactor. The two modes agree to
             well within a percent.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_radius_of_gyration(calculation_mode)



    # .....................................................................................
    #        
    def get_mean_end_to_end_distance(self, calculation_mode='scaling law'):
        """
        Returns the mean end-to-end distance (:math:`R_e`).

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
           Value equal to the average end-to-end distance (as defined by ``calculation_mode``).

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_end_to_end_distance(calculation_mode)

    # .....................................................................................
    #        
    def get_mean_hydrodynamic_radius(self, calculation_mode='kirkwood-riseman'):
        """
        Returns the average hydrodynamic radius, calculated either using the Kirkwood-Riseman
        equation or using the empirical Rg-to-Rh conversion scheme developed by Nygaard et al.

        In "kirkwood-riseman" mode the hydrodynamic radius is

        .. math::

           R_h = \\left\\langle \\frac{1}{r_{ij}} \\right\\rangle_{i \\neq j}^{-1}

        where the average runs over every pair of residues in the chain and, for each
        pair, over the AFRC's Gaussian inter-residue distance distribution. That inner
        average has the closed form :math:`\\langle 1/r_{ij} \\rangle = \\sqrt{6 / (\\pi
        \\langle r_{ij}^2 \\rangle)}`, so the result is exact for the model - no numerical
        integration is involved. Note that this is the mean of the *inverse* distance, as
        the Kirkwood-Riseman equation requires; it is not the inverse of the mean
        distance (for a Gaussian chain the two differ by a factor of :math:`4/\\pi`).
        This is the same form of the equation used by Nygaard et al. [1] and Pesce et al.
        [3], and by the ``mode='kr'`` option in SOURSOP.

        In "nygaard" mode the empirical relationship between :math:`R_g`, :math:`R_h` and
        chain length from Nygaard et al. [1] is applied to the mean radius of gyration.

        Parameters
        ----------

        calculation_mode : string (keyword)
            Defines how the hydrodynamic radius should be calculated. Must be one of either
            "kirkwood-riseman" or "nygaard".

        Returns
        -------
        float
           Value equal to the average hydrodynamic radius (in Angstroms).

        Raises
        ------
        AFRCException
           If the chain has fewer than two residues. In "kirkwood-riseman" mode there
           are then no inter-residue distances to average over, and in "nygaard" mode
           the :math:`N^{0.60} - N^{0.33}` denominator of the empirical relationship
           vanishes.

        References
        -------------
        [1] Nygaard M, Kragelund BB, Papaleo E, Lindorff-Larsen K. An Efficient
        Method for Estimating the Hydrodynamic Radius of Disordered Protein
        Conformations. Biophys J. 2017;113: 550–557.

        [2] Kirkwood, J. G., & Riseman, J. (1948). The Intrinsic Viscosities
        and Diffusion Constants of Flexible Macromolecules in Solution.
        The Journal of Chemical Physics, 16(6), 565–573.

        [3] Pesce, F., Newcombe, E. A., Seiffert, P., Tranchant, E. E.,
        Olsen, J. G., Grace, C. R., Kragelund, B. B., & Lindorff-Larsen, K.
        (2023). Assessment of models for calculating the hydrodynamic radius
        of intrinsically disordered proteins. Biophysical Journal, 122(2),
        310-321.


        """

        calculation_mode = validate_keyword(['kirkwood-riseman','nygaard'], calculation_mode, 'calculation_mode')

        # neither estimator is defined for a single residue: Kirkwood-Riseman has no
        # residue pairs to average over, and the Nygaard denominator N^0.6 - N^0.33 is
        # zero at N = 1 (which previously returned -0.0 with a divide-by-zero warning)
        n = len(self)
        if n < 2:
            raise AFRCException('The hydrodynamic radius requires a chain of at least two residues')

        if calculation_mode == 'nygaard':

            alpha1 = 0.216
            alpha2 = 4.06
            alpha3 = 0.821

            # first compute the rg
            rg = self.get_mean_radius_of_gyration()

            # precompute
            N_033 = np.power(n, 0.33)
            N_060 = np.power(n, 0.60)

            Rg_over_Rh = ((alpha1*(rg - alpha2*N_033)) / (N_060 - N_033)) + alpha3

            return (1/Rg_over_Rh)*rg

        elif calculation_mode == 'kirkwood-riseman':

            # the Kirkwood-Riseman equation averages the INVERSE inter-residue
            # distance, <1/r_ij>, over every pair of residues. Each PolymerObject in
            # the inter-residue matrix knows its own <1/r> exactly (closed form for
            # the Gaussian distribution), so we simply average those over the
            # non-redundant pairs and invert. Note this deliberately does not use
            # the distance map: 1/<r_ij> is a different quantity, and for a Gaussian
            # chain over-estimates Rh by a factor of 4/pi
            self.__build_matrix()

            inverse_distances = []
            for i in range(0, n):
                for j in range(i+1, n):
                    inverse_distances.append(self.matrix[i][j].get_mean_inverse_end_to_end_distance())

            return 1.0/np.mean(inverse_distances)

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

        tuple of np.ndarray
           A 2-pair tuple ``(distances, probabilities)`` where the first array is the
           distance (in Angstroms) and the second is the corresponding probability.

        """

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)
 
        if R1 == R2:
            return (np.array([0.0]), np.array([1.0]))

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
             then the complete Re distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Re = R0 * N^{0.5} is used.


        Returns
        -------

        float
           The mean distance (in Angstroms) between residues R1 and R2.


        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        if R1 == R2:
            return 0.0


        self.__build_matrix()

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

        R1 : int
           The first residue of the pair being investigated.

        R2 : int
           The second residue of the pair being investigated.

        calculation_mode : string (keyword)
             calculation_mode defines the mode in which the average is calculated, and can be
             set to either 'scaling law' (default) or 'distribution'. If 'distribution' is used
             then the complete Rg distribution is used to calculate the expected value. If the
             'scaling law' is used then the standard Rg = RG_R0 * N^{0.5} is used.


        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        self.__build_matrix()

        if R1 == R2:
            return 0.0

        return self.matrix[R1][R2].get_mean_radius_of_gyration(calculation_mode)
        

    # .....................................................................................
    #
    def sample_radius_of_gyration_distribution(self,n=1000):
        """
        Subsamples from the :math:`R_g` distribution to generate an uncorrelated 'trajectory'
        of points. Useful for creating a size-matched sample to compare with simulation
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
    def sample_end_to_end_distribution(self,n=1000):
        """
        Subsamples from the end-to-end distance distribution to generate an uncorrelated 
        'trajectory' of points. Useful for creating a size-matched sample to compare with
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

    
    # .....................................................................................
    #
    def sample_inter_residue_distance_distribution(self, R1, R2, n=1000):
        """
        Subsamples from the inter-residue distance distribution (between residues
        R1 and R2) to generate an uncorrelated 'trajectory' of points. Useful for
        creating a size-matched sample to compare with simulation data.

        Parameters
        ----------
        R1 : int
           The first residue of the pair being investigated.

        R2 : int
           The second residue of the pair being investigated.

        n : int
           Number of random values to sample (default = 1000)

        Returns
        -------
        np.ndarray
           Returns an n-length array with n independent values (floats)

        Raises
        ------
        AFRCException
           If either residue index is invalid.

        """

        # validate the indices before they are used - without this a negative index
        # silently wraps round and samples the wrong pair of residues, and an
        # out-of-range index raises a bare IndexError
        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        # construct the internal matrix of polymers
        self.__build_matrix()

        return self.matrix[R1][R2].sample_end_to_end_distribution(dist_size=n)


    # .....................................................................................
    #
    def get_contact_fraction(self, R1, R2, threshold):
        """
        Function that - given two residues (R1, and R2) and a distance threshold in Angstroms (threshold)
        returns the fraction of the time the center-of-mass distance between R1 and R2 is < threshold.

        Practically, if we set threshold = 5, this gives you the expected contact fraction for two residues,
        which is a useful normalization factor.

        Parameters
        -------------
        R1 : int
            First residue - must be between 0 and length of the polymer

        R2 : int
            Second residues - must also be between 0 and length of the polymer

        threshold : float
            A distance threshold in angstroms - can be a float or an int

        Returns
        ----------
        float
            Returns a single value between 0 and 1 that reports on the fraction
            of the time residues R1 and R2 are closer than ``threshold`` angstroms
            apart.

        """

        # validate the indices up front so that an invalid pair is rejected even when
        # R1 == R2 (which otherwise short-circuits below)
        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        if R1 == R2:
            return 1.0

        # get the distribution of inter-residue distances
        (r, p) = self.get_interresidue_distance_distribution(R1, R2)

        # p is a normalized probability MASS function (each bin holds the weight
        # of the interval [r, r + dr)), so the contact fraction is simply the
        # cumulative weight in the bins that lie below the threshold. Note this
        # was previously evaluated with the trapezoid rule, which halves the two
        # end bins and so systematically under-counted by half a bin's weight -
        # a ~1.5% relative error at a 5 A threshold on the default 0.05 A grid,
        # and considerably worse on the coarser adaptable grid
        return float(np.sum(p[r < threshold]))
        
    # .....................................................................................
    #
    def get_contact_map(self, threshold, symmetric_map=False):
        """
        Function that returns a contact map for the protein, where the contact map is a 
        square matrix where each element is the contact fraction between two residues. 

        Parameters
        -------------
        threshold : float
            A distance threshold in angstroms - can be a float or an int.

        symmetric_map : bool (default = False)
            If True, a full [n x n] matrix is returned, if False only the upper right triangle
            is returned.

        Returns
        ----------
        np.ndarray
            Returns a square matrix where each element is the contact fraction between 
            two residues. 

        """

        # initialize the contact map
        contact_map = np.zeros((len(self.seq),len(self.seq)))

        # for each pair of residues calculate the contact fraction
        for i in range(0,len(self.seq)):
            for j in range(i,len(self.seq)):
                contact_map[i,j] = self.get_contact_fraction(i, j, threshold)
                if symmetric_map:
                    contact_map[j,i] = contact_map[i,j]

        return contact_map
        
        


        
    def get_pre_profile(self, label_position, tau_c=4, t_delay=12, R_2D=14, W_H=2*np.pi*600e6, sample_size=10000):
        """
        Calculate the hypothetical paramagnetic relaxation enhancement (PRE) profile
        expected if a spin label were placed at position label_position. The
        only required input is the label position, but additional experimental
        parameters can be passed in as well.

        It's important to remember this method does not consider the explicit
        position of a spin label linker, but does provide a reference model
        as to the expected PRE profile if the chain behaved as an AFRC chain.
        

        Parameters
        -----------------------
        label_position : int
            Position along the chain that is labelled. Must be between 0 and the
            length of the sequence minus one.

        tau_c : float
            tau_c is the effective correlation time, measured in nanoseconds, 
            which is typically between 1 and 30. Default = 4

        t_delay : float
            Total duration of the INEPT delays from the PRE experiment, as 
            measured in ms. This will depend on the pulse sequence used, 
            but is typically around 1-30 ms for HSQC. Default = 12

        R_2D : float
            Is the transverse relaxation rate of the backbone amide protons in
            the diamagnetic form of the protein, measured in Hertz (i.e. 'per
            second'). A value of around 10 might be expected. Default = 14

        W_H : float
            Is the proton Larmor frequency as an *angular* frequency, in rad/s
            (:math:`\\omega_H = 2\\pi\\nu_H`). The spectral-density term of the
            Solomon-Bloembergen prefactor evaluates
            :math:`3\\tau_c/(1 + \\omega_H^2 \\tau_c^2)`, so the value passed here
            is used directly as the angular frequency. For a 600 MHz magnet pass
            ``2*np.pi*600e6`` (approximately 3.77e9), **not** ``600000000``;
            passing the linear frequency makes the dispersive term a factor of
            :math:`(2\\pi)^2 \\approx 39.5` too small, over-estimating gamma by
            roughly 10% at the default tau_c. This is the same convention used by
            SOURSOP's ``SSPRE`` class. Default = 2*pi*600e6 (a 600 MHz magnet)

        sample_size : int
            Number of conformations drawn from each inter-residue distance
            distribution when computing the relaxation rates. Larger values give a
            smoother profile at the cost of more compute. Default = 10000

        Returns
        -----------------------
        list
            Returns a 3-element list.

            [0] -  residue indices (starting at 0)
            [1] -  PRE profile (a value between 0 and 1)
            [2] -  PRE H1 relaxation profile (gamma), one array of per-conformation relaxation rates per residue

        Raises
        -----------------------
        AFRCException
            If label_position is not a valid residue index.

        References
        -------------
        [1] Meng, W., Lyle, N., Luan, B., Raleigh, D.P., and Pappu, R.V. 
        (2013). Experiments and simulations show how long-range
        contacts can form in expanded unfolded proteins with negligible 
        secondary structure.
        Proc. Natl. Acad. Sci. U. S. A. 110, 2123-2128.

        [2] Das, R.K., Huang, Y., Phillips, A.H., Kriwacki, R.W., and Pappu, 
        R.V. (2016). Cryptic sequence features within the disordered protein 
        p27Kip1 regulate cell cycle signaling. Proc. Natl. Acad. Sci. U. S. A. 
        113, 5616- 5621.
        
        [3] Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel, 
        O., & Raleigh, D. P. (2019). Unfolded states under folding conditions 
        accommodate sequence-specific conformational preferences with random 
        coil-like dimensions. Proceedings of the National Academy of Sciences 
        of the United States of America, 116(25), 12301–12310.

        [4] Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V., & 
        Holehouse, A. S. (2023). SOURSOP: A Python package for the analysis 
        of simulations of intrinsically disordered proteins. bioRxiv : The 
        Preprint Server for Biology. https://doi.org/10.1101/2023.02.16.528879


        """

        # this raises an AFRCException for a negative, out-of-range or non-integer
        # label position (rather than a TypeError for e.g. a string)
        label_position = self.__validate_residue_index(label_position)
        
        # local constants (show in a couple of units for clarity...)
        original_K = 1.2300e-32       # K constant in cm6*s-2
        K_IN_NM6   = original_K*1e42  # K constant in nm6 s-2

        # # convert tau_c to seconds and calculate tau_c squared
        tau_c = float(tau_c)/1000000000     # tau c in seconds
        tau_c_squared = tau_c * tau_c       #

        t_delay_in_seconds = t_delay/1000

        # compute the prefactor term which will be used when computing the PRE dependent 
        # relaxation profiles
        W_H_SQUARED = W_H*W_H
        PREFACTOR = (3 * tau_c)/(1 + W_H_SQUARED * tau_c_squared)
        PREFACTOR = (4*tau_c + PREFACTOR)
        PREFACTOR = PREFACTOR * K_IN_NM6

        gamma = []

        # finally for each residue calculate the r^6 distances associated with each frame and for EACH FRAME calculate
        # the PREFACTOR / r^6 value and then take the mean. THIS gives a different answer to if you take the mean distance
        # and calculate the PREFACTOR/<R^6> value because there is a non-linear mapping between relaxation and distance so
        # it's important the former method is used (i.e. only average at the end). This calculates the gamma coefficient for
        # each residue, which measures relaxation
        for idx in range(0, len(self.seq)):

            distances_nm = 0.1*self.sample_inter_residue_distance_distribution(label_position, idx, n=sample_size)
            distances_nm_r6 = np.power(distances_nm, 6)

            # old method - average here
            #gamma.append(np.mean(PREFACTOR/distances_nm_r6))

            # new method - take the full distribution of corrected values. Note that
            # for idx == label_position the distances are all zero (a residue with
            # itself), giving an infinite relaxation rate; this is expected and the
            # divide-by-zero is suppressed here.
            with np.errstate(divide='ignore'):
                gamma.append(PREFACTOR/distances_nm_r6)

        profile = []
        for g in gamma:

            # old method
            #profile.append((R_2D * np.exp(-g*t_delay_in_seconds)) / (R_2D + g))

            # new method
            profile.append(np.mean((R_2D * np.exp(-g*t_delay_in_seconds)) / (R_2D + g)))

        indices = np.arange(0,len(self.seq))
        return [indices, profile, gamma]
            


        

            
            
            




