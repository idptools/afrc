import numpy as np
from afrc.config import P_OF_R_RESOLUTION

class SAWException(Exception):
    pass

class SAW:
    """
    This class generates an object that returns polymer statistics consistent with a
    self-avoiding random walk (SAW). This model was developed by Jhullian 'J' Alston,
    and is based on the reference implementation by O'Brien et al [1]

    The chain is fixed at the good-solvent scaling exponent, held on the object as
    ``self.nu = 0.598``. That single value sets both the chain-length dependence of the
    end-to-end distance (:math:`R_{ee} = \\texttt{prefactor}\\,N^{\\nu}`) and the universal
    :math:`R_g/R_e` ratio, so the two are guaranteed to describe the same chain. To vary
    the exponent, use :class:`~afrc.polymer_models.nudep_saw.NuDepSAW` instead.

    [1] O’Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009).
    How accurate are polymer models in the analysis of Forster resonance
    energy transfer experiments on proteins? The Journal of Chemical Physics,
    130(12), 124903.

    [2] Le Guillou, J. C., & Zinn-Justin, J. (1977). Critical Exponents for the n-Vector
    Model in Three Dimensions from Field Theory. Physical Review Letters, 39(2), 95-98.

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

        # normalization constants for the des Cloizeaux scaling form, as tabulated by
        # O'Brien et al. These are a matched set: they are exactly the values required
        # to normalize P(r) and to set its root-mean-square to Ree given theta = 0.3
        # and delta = 2.5, so they should only ever be changed together.
        self.a = 3.67853
        self.b = 1.23152

        self.theta = 0.3
        self.delta = 2.5

        # the Flory scaling exponent and the (Le Guillou / Zinn-Justin) gamma exponent
        # for a self-avoiding walk. nu is used BOTH to set the chain-length dependence of
        # Ree (see __compute_end_to_end_distribution) and in the universal Rg/Re ratio
        # (see get_mean_radius_of_gyration) - it is held here as a single value so that
        # the two cannot drift apart
        self.nu = 0.598
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

        The mean is computed by integrating over the :math:`P(r)` vs. :math:`r`
        distribution (i.e. :math:`\\sum r \\cdot P(r)`), consistent with the
        convention used by the other models in this package.

        By default this uses a prefactor of 5.5 A (0.55 nanometers).

        Parameters
        ----------
        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the mean end-to-end distance.

        """

        [a, b] = self.get_end_to_end_distribution(prefactor)

        return np.sum(a * b)

    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self, prefactor=5.5):
        """
        Returns the root-mean-square end-to-end distance (:math:`\\sqrt{\\langle R_e^2 \\rangle}`)
        as calculated from the SAW model.

        The value is computed by taking the square root after integrating over
        :math:`P(r)` vs. :math:`r^2`.

        Parameters
        ----------
        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the root-mean-square end-to-end distance.

        """

        [a, b] = self.get_end_to_end_distribution(prefactor)

        return np.sqrt(np.sum(np.power(a, 2) * b))

    # .....................................................................................
    #        
    def get_mean_radius_of_gyration(self, prefactor=5.5):
        """
        Returns the mean radius of gyration (:math:`R_g`) for the SAW model.

        :math:`R_g` is obtained from the mean-squared end-to-end distance via the
        analytical ratio

        .. math::

           \\frac{\\langle R_g^2 \\rangle}{\\langle R_e^2 \\rangle} =
              \\frac{\\gamma(\\gamma + 1)}{2(\\gamma + 2\\nu)(\\gamma + 2\\nu + 1)}

        expressed in terms of the gamma exponent and the scaling exponent
        :math:`\\nu` (see [1]). Both are taken from the object's ``gamma`` and ``nu``
        attributes, so the exponent used here is by construction the same one that sets
        the chain-length dependence of :math:`R_{ee}`.

        Parameters
        ----------
        prefactor : float
            Prefactor that tunes the SAW dimensions. Default is 5.5 A.

        Returns
        -------
        float
           Value equal to the mean radius of gyration.

        """
        gamma = self.gamma
        nu = self.nu
        top = gamma*(gamma + 1)
        bottom = 2*(gamma + 2*nu)*(gamma + 2*nu + 1)

        # the ratio above relates mean-squared radii, so we use the
        # root-mean-square end-to-end distance (sqrt(<Re^2>)) here
        Ree = self.get_root_mean_squared_end_to_end_distance(prefactor=prefactor)

        return np.sqrt(Ree**2*(top/bottom))


    
    # .....................................................................................
    #        
    def __compute_end_to_end_distribution(self, prefactor):
        """
        Defines the end-to-end distribution based on the SAW as defined by
        https://aip.scitation.org/doi/10.1063/1.3082151
        . This is where we actually perform the polymer model calculation.


        """

        # a zero-length chain has all its weight at r = 0
        if self.zero_length:
            self.__p_of_Re_R = np.array([0.0])
            self.__p_of_Re_P = np.array([1.0])
            return

        # define the chainlength-dependent prefactor
        Ree = prefactor*np.power(self.nres, self.nu)

        # r values on a P(r) vs. r plot. We use the same grid as the parent AFRC model,
        # but ensure it always extends to at least 4*Ree - because Ree scales as N^nu
        # while the AFRC grid scales as N^0.5, for long chains the grid would otherwise
        # cut into the tail of the distribution and bias the mean and RMS values low
        upper = max(3*(7*np.power(self.nres, 0.5)), 4*Ree)
        p_dist = np.arange(0, upper, self.p_of_r_resolution)

        # define SAW normalization factors as defined by
        # https://aip.scitation.org/doi/10.1063/1.3082151
        a = self.a
        b = self.b

        theta = self.theta
        delta = self.delta

        # compute p(r) across the whole grid
        P_r_one = a/Ree
        P_r_two = np.power(p_dist/Ree, theta+2)
        P_r_three = np.exp(-b*np.power(p_dist/Ree, delta))

        p_val_raw = P_r_one*P_r_two*P_r_three

        # finally normalize so sums to 1.0 and assign to the object
        self.__p_of_Re_P = p_val_raw/np.sum(p_val_raw)
        self.__p_of_Re_R = p_dist
