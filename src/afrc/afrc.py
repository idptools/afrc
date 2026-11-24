"""
afrc.py

The Analytical Flory Random Coil (AFRC): a closed-form, sequence-specific
description of a polypeptide behaving as an ideal chain. It is fit to numerical
Flory Random Coil ensembles generated with the rotational isomeric state
approximation of Flory and Volkenstein, using backbone dihedral maps.

This module provides ``AnalyticalFRC``, the main user-facing object.

Copyright Alex Holehouse 2018-2026 (holehouselab.com).

"""
import numpy as np
from numpy.typing import NDArray

from .polymer import PolymerObject
from .config import P_OF_R_RESOLUTION, AA_list, RIJ_RMS_R0
from .ensemble import gaussian_chain_factor, sample_gaussian_chain, save_conformations, validate_n_conformations
from .ensemble_report import EnsembleReport, compare_ensemble_to_gaussian_model
from .exceptions import AFRCException
from .iofunctions import validate_keyword
from .polymer_models import wlc

# AFRCException is defined in exceptions.py but re-exported here so that the
# long-standing `from afrc.afrc import AFRCException` import path keeps working
__all__ = ['AFRCException', 'AnalyticalFRC']


class AnalyticalFRC:
    """
    The Analytical Flory Random Coil for a single amino acid sequence.

    This is the main object the package provides. You build it from one
    sequence and then ask it polymer questions - mean dimensions, full
    distributions, inter-residue distances, contact fractions, PRE profiles and
    so on.

    .. code-block:: python

       from afrc import AnalyticalFRC

       protein = AnalyticalFRC('KFGGPRDQGSRHDSEQDNSDNNTIFVQGLG')
       protein.get_mean_radius_of_gyration()

    Attributes
    ----------
    seq : str
        The sequence (upper-case).

    p_of_r_resolution : float
        Grid spacing (in Angstroms) used for every distribution.

    worm_like_chain : WormLikeChain
        A Zhou worm-like chain with the same number of residues (default
        ``lp`` and ``aa_size``), provided for convenience.

    Notes
    -----
    Nothing is computed until it is requested, so building the object is
    cheap. Methods that visit every pair of residues (distance maps, contact
    maps, internal scaling and the Kirkwood-Riseman :math:`R_h`) build an
    [n x n] matrix of sub-chains on first use, which is the most expensive
    step, and then keep it. Methods that need one pair, or one row of pairs
    (the PRE profile), build only the segments they need.

    Inter-residue quantities between residues :math:`i` and :math:`j` treat the
    segment between them as a chain of :math:`|i - j|` residues, taking its
    composition from the sequence between the two. Whole-chain quantities use
    all :math:`N` residues, so ``get_mean_interresidue_distance(0, N-1)`` is
    slightly smaller than ``get_mean_end_to_end_distance()``.

    """


    # .....................................................................................
    #
    def __init__(self, seq: str, adaptable_P_res: bool = False) -> None:
        """
        Create an AnalyticalFRC object from an amino acid sequence.

        Parameters
        ----------
        seq : str
            Amino acid sequence (case insensitive). Only the 20 standard amino
            acids are allowed.

        adaptable_P_res : bool
            If False (default) every distribution uses a fixed grid spacing of
            0.05 A. If True the spacing is set to :math:`d_{max}/500`, where
            :math:`d_{max} = 3.7N` approximates the contour length. This changes
            the discretization only, not the model.

        Raises
        ------
        AFRCException
            If ``seq`` is not a string, is empty, or contains a non-standard
            amino acid.

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

        # covariance factor for generating 3D conformations; built on first use by
        # sample_conformations()
        self.__ensemble_factor = None

        # finally we define other polymer models which are attached as their own
        # class objects
        self.worm_like_chain = wlc.WormLikeChain(seq, self.p_of_r_resolution)




    # .....................................................................................
    #
    def __check_seq_is_valid(self, seq):
        """
        Check that a sequence contains only the 20 standard amino acids.

        Parameters
        ----------
        seq : str
            Upper-case amino acid sequence.

        Raises
        ------
        AFRCException
            If any character is not a standard amino acid.

        """
        for i in seq:
            if i not in AA_list:
                raise AFRCException(f'Passed amino acid sequence contains non-standard amino acids [{i}]')



    # .....................................................................................
    #
    def __build_matrix(self):
        """
        Build the [n x n] matrix of inter-residue PolymerObjects, if not already built.

        Element ``[i][j]`` (and ``[j][i]``) describes the segment between
        residues ``i`` and ``j``, built from ``seq[i:j]``. The diagonal holds
        zero-length PolymerObjects. This is only needed by the methods that
        visit every pair (distance and contact maps, internal scaling and the
        Kirkwood-Riseman :math:`R_h`), so it is built once, on first use, and
        then kept.

        Note that the matrix entries do *not* cache their distributions. Keeping
        a full P(r) grid for every pair costs gigabytes of memory for a few
        hundred residues (6.8 GB at 400 residues), whereas recomputing one takes
        tens of microseconds.

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
                    self.matrix[i][j] = PolymerObject(subseq, self.p_of_r_resolution, cache_distributions=False)
                    self.matrix[j][i] = self.matrix[i][j]
        else:
            pass



    # .....................................................................................
    #
    def __get_pair(self, R1, R2):
        """
        Return the PolymerObject for the segment between two residues.

        If the inter-residue matrix has already been built we use its entry.
        Otherwise we build just this one segment - building the whole matrix
        for a single pair made one-off queries on long sequences very slow (over
        10 s for 1000 residues).

        Parameters
        ----------
        R1 : int
            Index of the first residue (already validated).

        R2 : int
            Index of the second residue (already validated). The order of ``R1``
            and ``R2`` does not matter.

        Returns
        -------
        PolymerObject
            The segment between the two residues, built from
            ``seq[min(R1, R2):max(R1, R2)]``, exactly as in the matrix. It does
            not cache its distributions.

        """

        if self.matrix is not False:
            return self.matrix[R1][R2]

        return PolymerObject(self.seq[min(R1, R2):max(R1, R2)], self.p_of_r_resolution, cache_distributions=False)



    # .....................................................................................
    #
    def __len__(self):
        """
        Return the number of residues in the sequence.

        Returns
        -------
        int
            Sequence length.

        """
        return len(self.seq)


    def __validate_residue_index(self, R):
        """
        Check that a residue index is a valid, zero-based position in the sequence.

        Parameters
        ----------
        R : int
            Residue index (or anything that can be cast to an int).

        Returns
        -------
        int
            The index, cast to an int.

        Raises
        ------
        AFRCException
            If ``R`` cannot be cast to an int, or lies outside
            ``[0, len(seq) - 1]``.

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
        Return the mean inter-residue distance for every pair of residues.

        Parameters
        ----------
        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle r_{ij} \\rangle = R_0 |i - j|^{0.5}`, or
            ``'distribution'``, which takes the expectation over each pair's
            distance distribution. The two agree to within about 0.3%.

        symmetric_map : bool
            If True, return the full symmetric matrix. If False (default), only
            the upper triangle is filled and the lower triangle is zero.

        Returns
        -------
        np.ndarray
            An [n x n] matrix of mean inter-residue distances (in Angstroms).

        Raises
        ------
        AFRCException
            If ``calculation_mode`` is not recognized.

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
        Return the internal scaling profile.

        This is the mean distance between residues that are :math:`|i - j|`
        apart in sequence, averaged over every such pair. A log-log fit of the
        profile gives a slope of 0.5.

        Parameters
        ----------
        calculation_mode : str
            Either ``'scaling law'`` (default) or ``'distribution'`` - see
            ``get_distance_map()``.

        Returns
        -------
        np.ndarray
            An [n-1 x 2] matrix with one row per sequence separation. The first
            column is :math:`|i - j|` (1 to n-1) and the second is the mean
            distance (in Angstroms) at that separation.

        Raises
        ------
        AFRCException
            If ``calculation_mode`` is not recognized.

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
        Return the radius of gyration (:math:`R_g`) distribution.

        This uses equation 3 of [Lhuillier1988]_ with the composition-weighted
        AFRC prefactor.

        Returns
        -------
        tuple of np.ndarray
            ``(radii, probabilities)``, where radii are in Angstroms and the
            probabilities are a normalized probability mass function (they sum
            to 1).

        """

        return self.full_seq_PO.get_radius_of_gyration_distribution()



    # .....................................................................................
    #
    def get_end_to_end_distribution(self):
        r"""
        Return the end-to-end distance (:math:`R_e`) distribution.

        This is the Gaussian chain distribution (see [Rubinstein2003]_)

        :math:`P(r) = 4\pi r^2 \left( \frac{3}{2\pi \langle r^2 \rangle} \right)^{3/2} e^{-\frac{3 r^2}{2 \langle r^2 \rangle}}`

        with :math:`\sqrt{\langle r^2 \rangle} = R_0^{rms} N^{0.5}`.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1).

        """

        return self.full_seq_PO.get_end_to_end_distribution()



    # .....................................................................................
    #
    def get_mean_radius_of_gyration(self, calculation_mode: str = 'distribution') -> float:
        """
        Return the mean radius of gyration, :math:`\\langle R_g \\rangle`.

        Parameters
        ----------
        calculation_mode : str
            Either ``'distribution'`` (default), which takes the expectation over
            the :math:`R_g` distribution, or ``'scaling law'``, which uses
            :math:`\\langle R_g \\rangle = R_0^{g} N^{0.5}` with the
            composition-weighted prefactor. The two agree to well within 0.1%.

        Returns
        -------
        float
            The mean radius of gyration (in Angstroms).

        Raises
        ------
        AFRCException
            If ``calculation_mode`` is not recognized.

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_radius_of_gyration(calculation_mode)



    # .....................................................................................
    #
    def get_mean_end_to_end_distance(self, calculation_mode: str = 'scaling law') -> float:
        """
        Return the mean end-to-end distance, :math:`\\langle R_e \\rangle`.

        Parameters
        ----------
        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle R_e \\rangle = R_0 N^{0.5}`, or ``'distribution'``,
            which takes the expectation over the end-to-end distribution. The two
            agree to within about 0.3%.

        Returns
        -------
        float
            The mean end-to-end distance (in Angstroms).

        Raises
        ------
        AFRCException
            If ``calculation_mode`` is not recognized.

        """

        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        return self.full_seq_PO.get_mean_end_to_end_distance(calculation_mode)

    # .....................................................................................
    #
    def get_mean_hydrodynamic_radius(self, calculation_mode='kirkwood-riseman'):
        """
        Return the mean hydrodynamic radius, :math:`R_h`.

        In ``'kirkwood-riseman'`` mode (default) this is

        .. math::

           R_h = \\left\\langle \\frac{1}{r_{ij}} \\right\\rangle_{i \\neq j}^{-1}

        where the average runs over every pair of residues and, for each pair,
        over its Gaussian distance distribution. The inner average has the closed
        form :math:`\\langle 1/r_{ij} \\rangle = \\sqrt{6 / (\\pi \\langle r_{ij}^2 \\rangle)}`,
        so the result is exact for the model. Note that this is the mean of the
        *inverse* distance, which is what the Kirkwood-Riseman equation needs - not
        the inverse of the mean distance, which for a Gaussian chain is larger by
        a factor of :math:`4/\\pi`. This is the same form used by Nygaard et al.
        [1], Pesce et al. [3] and SOURSOP's ``mode='kr'``.

        In ``'nygaard'`` mode we instead apply the empirical :math:`R_g`-to-
        :math:`R_h` conversion of Nygaard et al. [1] to the mean radius of
        gyration.

        Parameters
        ----------
        calculation_mode : str
            Either ``'kirkwood-riseman'`` (default) or ``'nygaard'``.

        Returns
        -------
        float
            The mean hydrodynamic radius (in Angstroms).

        Raises
        ------
        AFRCException
            If ``calculation_mode`` is not recognized, or the chain has fewer
            than two residues (Kirkwood-Riseman then has no pairs to average
            over, and the Nygaard denominator :math:`N^{0.60} - N^{0.33}` is
            zero).

        References
        ----------
        [1] Nygaard, M., Kragelund, B. B., Papaleo, E., & Lindorff-Larsen, K.
        (2017). An efficient method for estimating the hydrodynamic radius of
        disordered protein conformations. Biophysical Journal, 113(3), 550-557.

        [2] Kirkwood, J. G., & Riseman, J. (1948). The intrinsic viscosities
        and diffusion constants of flexible macromolecules in solution.
        The Journal of Chemical Physics, 16(6), 565-573.

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

        # calculation_mode has already been validated, so the only other option
        # is 'kirkwood-riseman'
        else:

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
        Return the distance distribution between two residues.

        Parameters
        ----------
        R1 : int
            Index of the first residue (zero-based).

        R2 : int
            Index of the second residue (zero-based). The order of ``R1`` and
            ``R2`` does not matter.

        Returns
        -------
        tuple of np.ndarray
            ``(distances, probabilities)``, where distances are in Angstroms and
            the probabilities are a normalized probability mass function (they
            sum to 1). If ``R1 == R2`` all the weight sits at zero.

        Raises
        ------
        AFRCException
            If either residue index is invalid.

        """

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        if R1 == R2:
            return (np.array([0.0]), np.array([1.0]))

        return self.__get_pair(R1, R2).get_end_to_end_distribution()



    # .....................................................................................
    #
    def get_mean_interresidue_distance(self, R1, R2, calculation_mode='scaling law'):
        """
        Return the mean distance between two residues.

        Parameters
        ----------
        R1 : int
            Index of the first residue (zero-based).

        R2 : int
            Index of the second residue (zero-based). The order of ``R1`` and
            ``R2`` does not matter.

        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle r_{ij} \\rangle = R_0 |i - j|^{0.5}`, or
            ``'distribution'``, which takes the expectation over the pair's
            distance distribution.

        Returns
        -------
        float
            The mean distance between ``R1`` and ``R2`` (in Angstroms), or 0.0 if
            ``R1 == R2``.

        Raises
        ------
        AFRCException
            If either residue index or ``calculation_mode`` is invalid.

        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        if R1 == R2:
            return 0.0

        return self.__get_pair(R1, R2).get_mean_end_to_end_distance(calculation_mode)


    # .....................................................................................
    #

    def get_mean_interresidue_radius_of_gyration(self, R1, R2, calculation_mode='scaling law'):
        """
        Return the mean radius of gyration of the segment between two residues.

        Parameters
        ----------
        R1 : int
            Index of the first residue (zero-based).

        R2 : int
            Index of the second residue (zero-based). The order of ``R1`` and
            ``R2`` does not matter.

        calculation_mode : str
            Either ``'scaling law'`` (default), which uses
            :math:`\\langle R_g \\rangle = R_0^{g} |i - j|^{0.5}`, or
            ``'distribution'``, which takes the expectation over the segment's
            :math:`R_g` distribution.

        Returns
        -------
        float
            The mean radius of gyration of the segment (in Angstroms), or 0.0 if
            ``R1 == R2``.

        Raises
        ------
        AFRCException
            If either residue index or ``calculation_mode`` is invalid.

        """
        calculation_mode = validate_keyword(['distribution','scaling law'], calculation_mode, 'calculation_mode')

        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        if R1 == R2:
            return 0.0

        return self.__get_pair(R1, R2).get_mean_radius_of_gyration(calculation_mode)


    # .....................................................................................
    #
    def sample_radius_of_gyration_distribution(self,n=1000):
        """
        Draw random radii of gyration from the :math:`R_g` distribution.

        Useful for building a size-matched, uncorrelated sample to compare
        against simulation data.

        Parameters
        ----------
        n : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``n`` independent radii of gyration (in Angstroms).

        """

        return self.full_seq_PO.sample_radius_of_gyration_distribution(dist_size=n)



    # .....................................................................................
    #
    def sample_end_to_end_distribution(self,n=1000):
        """
        Draw random end-to-end distances from the end-to-end distribution.

        Useful for building a size-matched, uncorrelated sample to compare
        against simulation data.

        Parameters
        ----------
        n : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``n`` independent end-to-end distances (in Angstroms).

        """

        return self.full_seq_PO.sample_end_to_end_distribution(dist_size=n)


    # .....................................................................................
    #
    def sample_inter_residue_distance_distribution(self, R1, R2, n=1000):
        """
        Draw random distances between two residues from their distance distribution.

        Useful for building a size-matched, uncorrelated sample to compare
        against simulation data.

        Parameters
        ----------
        R1 : int
            Index of the first residue (zero-based).

        R2 : int
            Index of the second residue (zero-based). The order of ``R1`` and
            ``R2`` does not matter.

        n : int
            Number of values to draw. Default is 1000.

        Returns
        -------
        np.ndarray
            ``n`` independent distances (in Angstroms). If ``R1 == R2`` every
            value is 0.

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

        return self.__get_pair(R1, R2).sample_end_to_end_distribution(dist_size=n)


    # .....................................................................................
    #
    def get_contact_fraction(self, R1, R2, threshold):
        """
        Return the fraction of the time two residues are closer than a threshold.

        This is the cumulative probability of the inter-residue distance
        distribution below ``threshold``. With a threshold of around 5 A this
        gives the expected contact fraction for a pair of residues under the
        AFRC, which is a useful normalization factor.

        Parameters
        ----------
        R1 : int
            Index of the first residue (zero-based).

        R2 : int
            Index of the second residue (zero-based).

        threshold : float
            Distance threshold (in Angstroms).

        Returns
        -------
        float
            The contact fraction, between 0 and 1. A residue is always in
            contact with itself, so ``R1 == R2`` returns 1.0.

        Raises
        ------
        AFRCException
            If either residue index is invalid or ``threshold`` is not a number.

        """

        # validate the indices up front so that an invalid pair is rejected even when
        # R1 == R2 (which otherwise short-circuits below)
        R1 = self.__validate_residue_index(R1)
        R2 = self.__validate_residue_index(R2)

        # the threshold must be a number; without this a string raised a bare
        # TypeError from inside NumPy
        try:
            threshold = float(threshold)
        except (TypeError, ValueError):
            raise AFRCException('Could not convert threshold [%s] to a number' % (threshold))

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
        Return the contact fraction for every pair of residues.

        Parameters
        ----------
        threshold : float
            Distance threshold (in Angstroms) - see ``get_contact_fraction()``.

        symmetric_map : bool
            If True, return the full symmetric matrix. If False (default), only
            the upper triangle (including the diagonal) is filled and the lower
            triangle is zero.

        Returns
        -------
        np.ndarray
            An [n x n] matrix of contact fractions, each between 0 and 1. The
            diagonal is 1.

        Raises
        ------
        AFRCException
            If ``threshold`` is not a number.

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
        Return the expected paramagnetic relaxation enhancement (PRE) profile for
        a spin label at ``label_position``.

        For every residue we draw ``sample_size`` distances to the label from the
        AFRC, convert each to a transverse relaxation rate

        .. math::

           \\Gamma_2 = \\frac{K}{r^6} \\left( 4\\tau_c + \\frac{3\\tau_c}{1 + \\omega_H^2 \\tau_c^2} \\right)

        with :math:`K = 1.23 \\times 10^{-32}` cm\\ :sup:`6` s\\ :sup:`-2`, and then
        average the intensity ratio
        :math:`R_{2D} e^{-\\Gamma_2 t} / (R_{2D} + \\Gamma_2)` over those
        samples. Averaging per conformer (rather than converting a mean distance)
        matters, because relaxation depends very non-linearly on distance.

        Note that the model does not account for the spin-label linker. It gives
        the profile expected if the chain behaved as an AFRC chain.

        Parameters
        ----------
        label_position : int
            Index of the labelled residue (zero-based).

        tau_c : float
            Effective correlation time, in nanoseconds. Typically between 1 and
            30. Default is 4.

        t_delay : float
            Total duration of the INEPT delays in the PRE experiment, in
            milliseconds. This depends on the pulse sequence, but is typically
            1-30 ms for an HSQC. Default is 12.

        R_2D : float
            Transverse relaxation rate of the backbone amide protons in the
            diamagnetic protein, in Hz (per second). Around 10 is typical.
            Default is 14.

        W_H : float
            Proton Larmor frequency as an *angular* frequency, in rad/s
            (:math:`\\omega_H = 2\\pi\\nu_H`). For a 600 MHz magnet pass
            ``2*np.pi*600e6`` (about 3.77e9), **not** ``600000000``; passing the
            linear frequency makes the dispersive term :math:`(2\\pi)^2 \\approx 39.5`
            times too small and over-estimates :math:`\\Gamma_2` by roughly 10% at
            the default ``tau_c``. This matches the convention used by SOURSOP's
            ``SSPRE`` class. Default is ``2*np.pi*600e6`` (a 600 MHz magnet).

        sample_size : int
            Number of distances drawn per residue. Larger values give a smoother
            profile at the cost of more compute. Default is 10000.

        Returns
        -------
        list
            A 3-element list:

            [0] - residue indices (np.ndarray, starting at 0)

            [1] - the PRE intensity-ratio profile (one value between 0 and 1 per
            residue). The labelled residue itself is 0.

            [2] - the per-conformer relaxation rates :math:`\\Gamma_2` (one
            np.ndarray of length ``sample_size`` per residue, in s\\ :sup:`-1`).
            These are infinite for the labelled residue.

        Raises
        ------
        AFRCException
            If ``label_position`` is not a valid residue index.

        References
        ----------
        [1] Meng, W., Lyle, N., Luan, B., Raleigh, D. P., & Pappu, R. V. (2013).
        Experiments and simulations show how long-range contacts can form in
        expanded unfolded proteins with negligible secondary structure.
        Proceedings of the National Academy of Sciences, 110(6), 2123-2128.

        [2] Das, R. K., Huang, Y., Phillips, A. H., Kriwacki, R. W., & Pappu,
        R. V. (2016). Cryptic sequence features within the disordered protein
        p27Kip1 regulate cell cycle signaling. Proceedings of the National
        Academy of Sciences, 113(20), 5616-5621.

        [3] Peran, I., Holehouse, A. S., Carrico, I. S., Pappu, R. V., Bilsel,
        O., & Raleigh, D. P. (2019). Unfolded states under folding conditions
        accommodate sequence-specific conformational preferences with random
        coil-like dimensions. Proceedings of the National Academy of Sciences,
        116(25), 12301-12310.

        [4] Lalmansingh, J. M., Keeley, A. T., Ruff, K. M., Pappu, R. V., &
        Holehouse, A. S. (2023). SOURSOP: A Python package for the analysis
        of simulations of intrinsically disordered proteins. Journal of
        Chemical Theory and Computation, 19(16), 5609-5620.

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


    # .....................................................................................
    #
    def get_mean_squared_distance_map(self) -> NDArray[np.float64]:
        """
        Return the mean-squared distance between every pair of residues.

        For residues :math:`i` and :math:`j` this is
        :math:`\\langle r_{ij}^2 \\rangle = (R_0^{rms})^2 |i - j|`, with
        :math:`R_0^{rms}` the average per-residue prefactor over the segment
        between them - exactly the value used by the inter-residue distance
        distributions. Because every AFRC inter-residue distance is Gaussian, this
        matrix fully specifies the model's pairwise statistics (and is what
        ``sample_conformations()`` draws from). Unlike ``get_distance_map()`` it
        does not build the inter-residue matrix; every pair is computed at once
        from a running sum of the per-residue prefactors, so it is fast even for
        long sequences.

        Returns
        -------
        np.ndarray
            Symmetric [n x n] matrix of mean-squared distances (in Angstroms
            squared), with zeros on the diagonal.

        """

        n = len(self.seq)
        prefactors = np.array([RIJ_RMS_R0[residue] for residue in self.seq])
        running_sum = np.concatenate(([0.0], np.cumsum(prefactors)))

        i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
        lower, upper = np.minimum(i, j), np.maximum(i, j)
        separation = upper - lower

        # (sum of prefactors)^2 / |i - j| = (mean prefactor)^2 * |i - j|; the
        # diagonal (zero separation) is zero
        segment_sum = running_sum[upper] - running_sum[lower]
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(separation > 0, segment_sum*segment_sum/separation, 0.0)


    # .....................................................................................
    #
    def sample_conformations(self, n: int = 1000, seed: int | np.random.Generator | None = None) -> NDArray[np.float64]:
        """
        Generate 3D conformations (one bead per residue) drawn from the AFRC.

        The AFRC defines every inter-residue distance as Gaussian, with a
        mean-squared value set by the composition of the sequence between the two
        residues. Taken together, those distances define a Gaussian chain whose
        bead coordinates are jointly Gaussian, with covariance

        .. math::

           C = -\\frac{1}{6} J D J, \\qquad D_{ij} = \\langle r_{ij}^2 \\rangle,
           \\qquad J = I - \\frac{1}{n}\\mathbf{1}\\mathbf{1}^T

        for each of x, y and z. We factor :math:`C` once per sequence and then
        draw conformations from it directly. Every inter-residue distance in the
        resulting ensemble has exactly the AFRC's distance distribution, for
        every pair of residues at once - this is not an approximation or a fit.

        There are a few things worth knowing about the ensembles:

        * Adjacent beads are not at a fixed 3.8 A spacing. Their separation
          follows the AFRC's own distribution for neighbouring residues (a mean of
          ~5.8 A, with most values between 2 and 10 A).
        * The AFRC is an ideal chain, so there is no excluded volume and beads can
          overlap.
        * The end-to-end distance of the ensemble is the distance between the
          first and last residues, which the AFRC treats as a chain of N-1
          residues. It therefore matches
          ``get_mean_interresidue_distance(0, N-1)``, which is slightly smaller
          than ``get_mean_end_to_end_distance()`` (N residues).
        * The radius of gyration is that of the Gaussian chain, which agrees with
          ``get_mean_radius_of_gyration()`` to within ~2% (the AFRC's Rg was
          calibrated separately from its distances).

        Each conformation is centred on the origin and randomly oriented.

        Parameters
        ----------
        n : int
            Number of conformations to generate. Default is 1000.

        seed : int, np.random.Generator or None
            Seed (or generator) for the random numbers, for reproducible
            ensembles. If None (default) a fresh, unpredictable seed is used.

        Returns
        -------
        np.ndarray
            Array of shape [n x N x 3] with the bead coordinates, in Angstroms,
            where N is the sequence length.

        Raises
        ------
        AFRCException
            If ``n`` is not a positive integer.

        Examples
        --------
        >>> protein = AnalyticalFRC('MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI')
        >>> xyz = protein.sample_conformations(n=5000, seed=1)
        >>> xyz.shape
        (5000, 50, 3)

        """

        n_conformations = validate_n_conformations(n)

        # factoring the covariance is O(N^3) but only needs doing once per sequence
        if self.__ensemble_factor is None:
            self.__ensemble_factor = gaussian_chain_factor(self.get_mean_squared_distance_map())

        return sample_gaussian_chain(self.__ensemble_factor, n_conformations, np.random.default_rng(seed))


    # .....................................................................................
    #
    def save_ensemble(self, filename: str, n: int = 1000, seed: int | np.random.Generator | None = None, pdb_only: bool = False) -> NDArray[np.float64]:
        """
        Generate an AFRC ensemble and write it to disk.

        This draws ``n`` conformations with ``sample_conformations()`` and, by
        default, writes them as a PDB/XTC pair: ``<filename>.pdb`` (the topology,
        one CA bead per residue, holding the first conformation) and
        ``<filename>.xtc`` (every conformation). Load them together with, for
        example, ``mdtraj.load('ens.xtc', top='ens.pdb')`` or SOURSOP's
        ``SSTrajectory('ens.xtc', 'ens.pdb')``. Writing the XTC file needs mdtraj
        (``pip install mdtraj``).

        With ``pdb_only=True`` every conformation is instead written to a single
        multi-model ``<filename>.pdb``, which needs no extra dependencies but is
        larger and limited to 9999 conformations.

        Parameters
        ----------
        filename : str
            Output path without an extension; ``.pdb`` (and ``.xtc``) are added.
            If the name already ends in ``.pdb`` or ``.xtc`` that extension is
            dropped first.

        n : int
            Number of conformations to generate. Default is 1000.

        seed : int, np.random.Generator or None
            Seed (or generator) for the random numbers, for reproducible
            ensembles. If None (default) a fresh, unpredictable seed is used.

        pdb_only : bool
            If True, write a single multi-model PDB file instead of a PDB/XTC
            pair. Default is False.

        Returns
        -------
        np.ndarray
            The conformations that were written, as an array of shape [n x N x 3]
            in Angstroms.

        Raises
        ------
        AFRCException
            If ``n`` is not a positive integer, mdtraj is not installed (for a
            PDB/XTC pair), the sequence is longer than the PDB format allows (9999
            residues), there are more than 9999 conformations (with
            ``pdb_only=True``), or a conformation that goes in the PDB file does
            not fit in the PDB coordinate format.

        """

        from . import __version__

        conformations = self.sample_conformations(n=n, seed=seed)
        save_conformations(conformations, self.seq, filename, pdb_only=pdb_only,
                           remark=f'AFRC ensemble from afrc {__version__}: {len(conformations)} conformations')
        return conformations


    # .....................................................................................
    #
    def check_ensemble(self, conformations: NDArray[np.float64]) -> EnsembleReport:
        """
        Check how well an ensemble reproduces the AFRC's statistics.

        Every quantity the AFRC fixes exactly is measured in the ensemble, with a
        standard error estimated from the ensemble, and compared with the model:
        the root-mean-square radius of gyration, the first-to-last bead distance,
        the Kirkwood-Riseman hydrodynamic radius, the mean-squared and mean
        inter-residue distances (in bands of sequence separation), and the full
        distribution of three representative distances. The AFRC's separately
        calibrated :math:`\\langle R_g \\rangle` and its whole-chain
        :math:`\\langle R_e \\rangle` are reported as context. This is the report
        that ``afrc-ensemble`` prints.

        Parameters
        ----------
        conformations : np.ndarray
            Array of shape [n_conformations x N x 3] (in Angstroms), e.g. from
            ``sample_conformations()``.

        Returns
        -------
        EnsembleReport
            The checks; ``report.passed`` says whether the ensemble is model-like
            and ``report.format()`` gives a printable report.

        Raises
        ------
        AFRCException
            If the conformations do not match this sequence.

        """
        return compare_ensemble_to_gaussian_model(conformations, self.get_mean_squared_distance_map(),
                                                  model_name='AFRC',
                                                  reference_mean_rg=self.get_mean_radius_of_gyration(),
                                                  reference_mean_re=self.get_mean_end_to_end_distance())
