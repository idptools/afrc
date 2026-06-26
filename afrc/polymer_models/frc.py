import numpy as np
from afrc.config import P_OF_R_RESOLUTION
from numpy.random import choice

class FRCException(Exception):
    pass

class FreelyRotatingChain:
    """
    This class generates an object that returns polymer statistics consistent with a
    freely rotating chain (FRC) of ``N`` bonds of length ``b`` with a fixed bond angle
    and unrestricted (free) torsion angles.

    The freely rotating chain is an ideal chain: like the Analytical Flory Random Coil it
    has Gaussian end-to-end statistics with a true scaling exponent of 0.5, but its size is
    set by a single stiffness parameter - the characteristic ratio :math:`C_\\infty`. The
    mean-squared end-to-end distance follows the exact finite-N freely-rotating-chain result

        ⟨R²⟩ = C∞·N·b² − 2 b² α (1 − α^N) / (1 − α)²,      α = (C∞ − 1)/(C∞ + 1),

    where α is the cosine of the angle between successive bonds. Setting ``c_inf = 1``
    (α = 0) recovers the freely jointed chain, while ``c_inf = 2`` corresponds to a
    tetrahedral backbone angle. Note that a freely rotating chain cannot reproduce the much
    larger characteristic ratio of a real polypeptide (:math:`C_\\infty \\approx 9`), which
    arises from hindered/restricted rotation - use the Analytical Flory Random Coil for that.

    This is a composition-independent model: the sequence is used only to set the number of
    bonds. It is included as an additional reference model.

    [1] Flory, P. J. (1969). Statistical Mechanics of Chain Molecules. Wiley-Interscience.

    [2] Rubinstein, M., & Colby, R. H. (2003). Polymer Physics. Oxford University Press.

    """

    # .....................................................................................
    #
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION, b=3.8, c_inf=2.0):
        """
        Method to create a FreelyRotatingChain object. Seq should be a valid upper-case
        amino acid sequence and p_of_r_resolution defines the resolution (in angstroms)
        to be used for distributions.

        By default p_of_r_resolution is taken from the config.py file in the afrc package
        which defines the resolution at 0.05 A.

        Parameters
        -----------
        seq : str
            Amino acid sequence (used only to calculate the number of bonds).

        p_of_r_resolution : float
            Bin width for building probability distributions. In Angstroms.

        b : float
            Bond (segment) length, in Angstroms. The default of 3.8 A corresponds to the
            Cα-Cα distance, i.e. one virtual bond per residue.

        c_inf : float
            Characteristic ratio :math:`C_\\infty`, a dimensionless measure of chain
            stiffness defined as :math:`(1 + \\alpha)/(1 - \\alpha)` where α is the cosine
            of the angle between successive bonds. ``c_inf = 1`` recovers the freely jointed
            chain; ``c_inf = 2`` corresponds to a tetrahedral bond angle. Must be > 0.

        """

        # cast and sanity check the free parameters
        self.b = float(b)
        if self.b <= 0:
            raise FRCException('Error, b (bond length) cannot be less than or equal to 0')

        self.c_inf = float(c_inf)
        if self.c_inf <= 0:
            raise FRCException('Error, c_inf (characteristic ratio) cannot be less than or equal to 0')

        # set sequence info - the number of bonds
        self.nres = len(seq)

        # p_of_r_resolution defines the P(r) resolution in angstroms - i.e. basically
        # the spacing between r values in a P(r) vs. r plot
        self.p_of_r_resolution = p_of_r_resolution

        # set distribution info to false - these are calculated if/when needed
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
        Defines the end-to-end distribution based on the freely rotating chain (FRC).

        The freely rotating chain is an ideal chain, so the distribution is Gaussian, but
        its size is set by the characteristic ratio. This is a composition independent
        model for which the end-to-end distance depends solely on the number of amino acids.
        It is included here as an additional reference model.

        Returns
        -------

        tuple of arrays
           A 2-pair tuple of numpy arrays where the first is the distance (in Angstroms) and
           the second array is the probability of that distance.

        """
        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self):
        """
        Returns the mean end-to-end distance (:math:`R_e`) for the freely rotating chain.

        The mean is computed by integrating over the :math:`P(r)` vs. :math:`r`
        distribution (i.e. :math:`\\sum r \\cdot P(r)`), consistent with the convention
        used by the other models in this package.

        Returns
        -------
        float
           Value equal to the mean end-to-end distance.

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sum(a * b)


    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self):
        """
        Returns the root-mean-square end-to-end distance
        (:math:`\\sqrt{\\langle R_e^2 \\rangle}`) for the freely rotating chain.

        The value is computed by taking the square root after integrating over
        :math:`P(r)` vs. :math:`r^2`. In the long-chain limit this approaches
        :math:`\\sqrt{C_\\infty}\\, b\\sqrt{N}`.

        Returns
        -------
        float
           Value equal to the root-mean-square end-to-end distance.

        """
        [a, b] = self.get_end_to_end_distribution()

        return np.sqrt(np.sum(b * np.power(a, 2)))


    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self):
        """
        Returns the mean radius of gyration (:math:`R_g`) for the freely rotating chain.

        For an ideal chain the radius of gyration is related to the root-mean-square
        end-to-end distance by :math:`R_g = \\sqrt{\\langle R_e^2 \\rangle / 6}`, and this
        relationship is used here.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        return self.get_root_mean_squared_end_to_end_distance() / np.sqrt(6)


    # .....................................................................................
    #
    def sample_end_to_end_distribution(self, n=1000):
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
        if self.zero_length:
            return np.repeat(0.0, n)

        if self.__p_of_Re_R is False:
            self.__compute_end_to_end_distribution()

        return choice(self.__p_of_Re_R, n, p=self.__p_of_Re_P)


    # .....................................................................................
    #
    def __compute_end_to_end_distribution(self):
        """
        Defines the end-to-end distribution based on the freely rotating chain (FRC).
        This is where we actually perform the polymer model calculation.

        The mean-squared end-to-end distance uses the exact finite-N freely-rotating-chain
        result, and the distribution is the corresponding Gaussian chain form

            P(r) = 4π r² (3 / 2π⟨R²⟩)^{3/2} exp( -3 r² / 2⟨R²⟩ ).

        """

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        N = self.nres
        b = self.b

        # the cosine of the angle between successive bonds, recovered from the
        # characteristic ratio: C∞ = (1 + α)/(1 - α)
        alpha = (self.c_inf - 1.0) / (self.c_inf + 1.0)

        # exact finite-N mean-squared end-to-end distance for the freely rotating chain.
        # The first term is the long-chain limit (C∞ N b²); the second is the finite-size
        # correction (which vanishes when α = 0, i.e. the freely jointed chain).
        mean_sq_re = self.c_inf * N * np.power(b, 2)
        if alpha != 0:
            mean_sq_re = mean_sq_re - 2 * np.power(b, 2) * alpha * (1 - np.power(alpha, N)) / np.power(1 - alpha, 2)

        # build an r-grid that comfortably captures the Gaussian tail regardless of the
        # chosen stiffness/segment length
        rms = np.sqrt(mean_sq_re)
        p_dist = np.arange(0, 4.0*rms, self.p_of_r_resolution)

        # standard Gaussian chain end-to-end distribution
        A = np.power(3.0/(2*np.pi*mean_sq_re), 1.5)
        p_val_raw = 4*np.pi*np.power(p_dist, 2)*A*np.exp(-(3*np.power(p_dist, 2))/(2*mean_sq_re))

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw / np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist
