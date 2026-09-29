"""
saw.py

Self-avoiding walk (SAW) model at a fixed good-solvent scaling exponent.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from afrc.config import P_OF_R_RESOLUTION

class SAWException(Exception):
    """Exception raised by the self-avoiding walk model."""
    pass

class SAW:
    """
    Self-avoiding walk, using the des Cloizeaux scaling form as implemented by
    O'Brien et al. [1]. This model was developed by Jhullian 'J' Alston.

    The scaling exponent is fixed at the good-solvent value, held on the object
    as ``self.nu = 0.598``. That single value sets both the chain-length
    dependence of the size scale (:math:`R_{ee} = \\texttt{prefactor}\\,N^{\\nu}`)
    and the universal :math:`R_g/R_e` ratio, so the two always describe the same
    chain. To vary the exponent, use
    :class:`~afrc.polymer_models.nudep_saw.NuDepSAW`.

    This is a composition-independent reference model: the sequence is only used
    to set the number of residues. The overall size is set by the ``prefactor``
    argument accepted by every method (default 5.5 A).

    References
    ----------
    [1] O'Brien, E. P., Morrison, G., Brooks, B. R., & Thirumalai, D. (2009).
    How accurate are polymer models in the analysis of Forster resonance
    energy transfer experiments on proteins? The Journal of Chemical Physics,
    130(12), 124903.

    [2] Le Guillou, J. C., & Zinn-Justin, J. (1977). Critical Exponents for the
    n-Vector Model in Three Dimensions from Field Theory. Physical Review
    Letters, 39(2), 95-98.

    """

    # .....................................................................................
    #
    def __init__(self, seq, p_of_r_resolution=P_OF_R_RESOLUTION):
        """
        Create a SAW object.

        Parameters
        ----------
        seq : str
            Amino acid sequence. Only its length is used, and it is not
            validated.

        p_of_r_resolution : float
            Grid spacing (in Angstroms) used for the distribution. Default is
            0.05 A.

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
        Return the end-to-end distance distribution.

        Because ``prefactor`` can change between calls, the distribution is
        recomputed every time rather than cached.

        Parameters
        ----------
        prefactor : float
            Size scale in Angstroms, such that the root-mean-square end-to-end
            distance is ``prefactor * N**0.598``. Default is 5.5 A; this should
            be tuned to match excluded-volume simulations for quantitative work.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1).

        Raises
        ------
        SAWException
            If ``prefactor`` is not positive.

        """

        self.__compute_end_to_end_distribution(prefactor)

        return (self.__p_of_Re_R, self.__p_of_Re_P)


    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self, prefactor=5.5):
        """
        Return the mean end-to-end distance, :math:`\\langle R_e \\rangle`.

        This is the expectation over the end-to-end distribution,
        :math:`\\sum r P(r)`.

        Parameters
        ----------
        prefactor : float
            Size scale in Angstroms (see ``get_end_to_end_distribution``).
            Default is 5.5 A.

        Returns
        -------
        float
            The mean end-to-end distance (in Angstroms).

        Raises
        ------
        SAWException
            If ``prefactor`` is not positive.

        """

        [a, b] = self.get_end_to_end_distribution(prefactor)

        return np.sum(a * b)

    # .....................................................................................
    #
    def get_root_mean_squared_end_to_end_distance(self, prefactor=5.5):
        """
        Return the root-mean-square end-to-end distance,
        :math:`\\sqrt{\\langle R_e^2 \\rangle}`.

        This is the square root of :math:`\\sum r^2 P(r)` over the end-to-end
        distribution, and equals ``prefactor * N**0.598``.

        Parameters
        ----------
        prefactor : float
            Size scale in Angstroms (see ``get_end_to_end_distribution``).
            Default is 5.5 A.

        Returns
        -------
        float
            The root-mean-square end-to-end distance (in Angstroms).

        Raises
        ------
        SAWException
            If ``prefactor`` is not positive.

        """

        [a, b] = self.get_end_to_end_distribution(prefactor)

        return np.sqrt(np.sum(np.power(a, 2) * b))

    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self, prefactor=5.5):
        """
        Return the root-mean-square radius of gyration,
        :math:`\\sqrt{\\langle R_g^2 \\rangle}`.

        This is obtained from the root-mean-square end-to-end distance via the
        universal ratio

        .. math::

           \\frac{\\langle R_g^2 \\rangle}{\\langle R_e^2 \\rangle} =
              \\frac{\\gamma(\\gamma + 1)}{2(\\gamma + 2\\nu)(\\gamma + 2\\nu + 1)}

        using the object's ``gamma`` and ``nu`` attributes, so the exponent here
        is by construction the one that sets the size scale. Note that despite
        the method name this is the *root-mean-square* radius of gyration, not
        :math:`\\langle R_g \\rangle`. The name is kept for consistency with the
        other models.

        Parameters
        ----------
        prefactor : float
            Size scale in Angstroms (see ``get_end_to_end_distribution``).
            Default is 5.5 A.

        Returns
        -------
        float
            The root-mean-square radius of gyration (in Angstroms).

        Raises
        ------
        SAWException
            If ``prefactor`` is not positive.

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
        Build the end-to-end distribution for a given prefactor.

        .. math::

           P(r) = \\frac{a}{R_{ee}} \\left( \\frac{r}{R_{ee}} \\right)^{2+\\theta}
                  \\exp\\left[ -b \\left( \\frac{r}{R_{ee}} \\right)^{\\delta} \\right],
           \\qquad R_{ee} = \\texttt{prefactor}\\,N^{\\nu}

        with :math:`\\theta = 0.3`, :math:`\\delta = 2.5`, :math:`a = 3.67853`
        and :math:`b = 1.23152`. The grid runs from 0 to the larger of
        :math:`21\\sqrt{N}` and :math:`4R_{ee}`.

        Parameters
        ----------
        prefactor : float
            Size scale in Angstroms.

        Raises
        ------
        SAWException
            If ``prefactor`` is not positive.

        """

        # the prefactor is a length scale, so it must be positive - a zero or negative
        # value previously produced an all-NaN distribution without complaint
        prefactor = float(prefactor)
        if prefactor <= 0:
            raise SAWException('Error, prefactor must be greater than 0 (it sets the size scale in Angstroms)')

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
