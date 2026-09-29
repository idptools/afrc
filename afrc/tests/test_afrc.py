"""
Tests for ``afrc/afrc.py`` - the user-facing ``AnalyticalFRC`` class.

Shared fixtures (``protein``, ``all_aa``, ``test_seq``) come from conftest.py.
"""

import sys

import numpy as np
import pytest

import afrc
from afrc.afrc import AFRCException


# ---------------------------------------------------------------------------
# import / construction / validation
# ---------------------------------------------------------------------------
def test_afrc_imported():
    """Sample test, will always pass so long as the import statement worked."""
    assert "afrc" in sys.modules


def test_construction_and_len(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    assert len(p) == len(all_aa)
    assert p.seq == all_aa


def test_construction_is_case_insensitive(all_aa):
    upper = afrc.AnalyticalFRC(all_aa)
    lower = afrc.AnalyticalFRC(all_aa.lower())
    assert lower.seq == all_aa
    assert lower.get_mean_radius_of_gyration() == upper.get_mean_radius_of_gyration()


def test_invalid_amino_acid_raises():
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC('ACDEFZ')


def test_non_string_input_raises():
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC(12345)


def test_adaptable_resolution(all_aa):
    p = afrc.AnalyticalFRC(all_aa, adaptable_P_res=True)
    expected = (3.7 * len(all_aa)) / 500.0
    assert p.p_of_r_resolution == pytest.approx(expected)


def test_default_resolution(all_aa):
    from afrc.config import P_OF_R_RESOLUTION
    p = afrc.AnalyticalFRC(all_aa)
    assert p.p_of_r_resolution == P_OF_R_RESOLUTION


# ---------------------------------------------------------------------------
# distributions are valid PMFs
# ---------------------------------------------------------------------------
def test_end_to_end_distribution_is_valid_pmf(protein):
    dist, prob = protein.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)
    assert dist[0] == pytest.approx(0.0)


def test_rg_distribution_is_valid_pmf(protein):
    dist, prob = protein.get_radius_of_gyration_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_interresidue_distribution_is_valid_pmf(protein):
    dist, prob = protein.get_interresidue_distance_distribution(10, 80)
    assert len(dist) == len(prob)
    assert np.sum(prob) == pytest.approx(1.0)


def test_interresidue_distribution_same_residue(protein):
    dist, prob = protein.get_interresidue_distance_distribution(5, 5)
    assert dist[0] == pytest.approx(0.0)
    assert prob[0] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# analytic / first-principles relationships
# ---------------------------------------------------------------------------
def test_internal_scaling_exponent(protein):
    """A Flory Random Coil should have an apparent scaling exponent of ~0.5."""
    iscaling = protein.get_internal_scaling()
    x = np.log(iscaling[1:, 0])
    y = np.log(iscaling[1:, 1])
    slope = np.polyfit(x, y, 1)[0]
    assert slope == pytest.approx(0.5, abs=0.01)


def test_internal_scaling_shape(protein):
    """One row per sequence separation |i-j| = 1 ... n-1, columns [|i-j|, <r>]."""
    iscaling = protein.get_internal_scaling()
    assert iscaling.shape == (len(protein) - 1, 2)
    assert np.array_equal(iscaling[:, 0], np.arange(1, len(protein)))


def test_rg_scaling_law_and_distribution_agree(protein):
    """The two ways of computing <Rg> should give effectively the same answer.

    The 'scaling law' mode uses the calibrated Rg prefactor (RG_R0), which is
    back-calculated from the same Lhuillier fits that generate the distribution, so
    the two routes agree to well within a percent.
    """
    rg_law = protein.get_mean_radius_of_gyration('scaling law')
    rg_dist = protein.get_mean_radius_of_gyration('distribution')
    assert rg_law == pytest.approx(rg_dist, rel=0.01)


def test_rg_is_close_to_re_over_sqrt6(protein):
    """Rg should sit near - but not exactly at - the ideal-chain Re/sqrt(6) value.

    Re/sqrt(6) is exact for *root-mean-square* radii; applied to the mean Re it
    under-estimates <Rg> by around 5%, which is why the scaling-law mode uses the
    calibrated prefactor instead.
    """
    re = protein.get_mean_end_to_end_distance('scaling law')
    rg = protein.get_mean_radius_of_gyration('scaling law')
    assert rg == pytest.approx(re / np.sqrt(6), rel=0.1)
    assert rg > re / np.sqrt(6)


def test_distance_distribution_and_scaling_law_agree(protein):
    """mean Re from the distribution and the scaling law should be close."""
    re_law = protein.get_mean_end_to_end_distance('scaling law')
    re_dist = protein.get_mean_end_to_end_distance('distribution')
    assert re_dist == pytest.approx(re_law, rel=0.05)


def test_invalid_calculation_mode_raises(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_end_to_end_distance('not-a-mode')


# ---------------------------------------------------------------------------
# residue-index validation (via public methods that use it)
# ---------------------------------------------------------------------------
def test_negative_residue_index_raises(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_distance(-1, 5)


def test_out_of_range_residue_index_raises(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_distance(0, len(protein))


def test_non_castable_residue_index_raises(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_distance('abc', 5)


def test_mean_interresidue_distance_same_residue(protein):
    assert protein.get_mean_interresidue_distance(7, 7) == 0.0


# ---------------------------------------------------------------------------
# distance & contact maps
# ---------------------------------------------------------------------------
def test_distance_map_shape_and_triangularity(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    dm = p.get_distance_map()
    n = len(all_aa)
    assert dm.shape == (n, n)
    # lower triangle should be zero when symmetric_map is False
    assert np.allclose(np.tril(dm, -1), 0.0)


def test_distance_map_symmetric(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    dm = p.get_distance_map(symmetric_map=True)
    assert np.allclose(dm, dm.T)


def test_contact_fraction_bounds_and_identity(protein):
    assert protein.get_contact_fraction(10, 10, 5) == 1.0
    frac = protein.get_contact_fraction(10, 40, 20.0)
    assert 0.0 <= frac <= 1.0


def test_contact_fraction_monotonic_in_threshold(protein):
    low = protein.get_contact_fraction(10, 40, 10.0)
    high = protein.get_contact_fraction(10, 40, 30.0)
    assert high >= low


def test_contact_fraction_below_grid_resolution(protein):
    """Regression: a threshold below one grid step raised an IndexError."""
    assert protein.get_contact_fraction(10, 40, 0.01) == 0.0
    assert protein.get_contact_fraction(10, 40, protein.p_of_r_resolution) == 0.0


def test_contact_fraction_validates_indices(protein):
    """Indices are checked even for the R1 == R2 short-circuit."""
    with pytest.raises(AFRCException):
        protein.get_contact_fraction(-1, -1, 5.0)
    with pytest.raises(AFRCException):
        protein.get_contact_fraction(0, len(protein), 5.0)


def test_contact_map_shape(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    cm = p.get_contact_map(15.0)
    n = len(all_aa)
    assert cm.shape == (n, n)
    assert np.all((cm >= 0) & (cm <= 1))


# ---------------------------------------------------------------------------
# hydrodynamic radius and PRE
# ---------------------------------------------------------------------------
def test_hydrodynamic_radius_modes(protein):
    rh_kr = protein.get_mean_hydrodynamic_radius('kirkwood-riseman')
    rh_ny = protein.get_mean_hydrodynamic_radius('nygaard')
    assert rh_kr > 0
    assert rh_ny > 0


def test_kirkwood_riseman_averages_inverse_distances(protein):
    """Regression: Rh was 1/<r_ij> (the inverse of the mean distance) rather than
    the Kirkwood-Riseman 1/<1/r_ij>; for a Gaussian chain those differ by 4/pi."""
    n = len(protein)
    dm = protein.get_distance_map()  # also builds the inter-residue matrix
    inv = []
    for i in range(n):
        for j in range(i + 1, n):
            rms = protein.matrix[i][j].RMS_Re_scaling
            inv.append(np.sqrt(6 / (np.pi * rms**2)))
    expected = 1 / np.mean(inv)
    rh = protein.get_mean_hydrodynamic_radius('kirkwood-riseman')
    assert rh == pytest.approx(expected, rel=1e-12)

    # the old estimator, for the record: inverse of the mean distance
    old = 1 / np.mean(1 / dm[dm != 0])
    assert old / rh == pytest.approx(4 / np.pi, rel=0.01)


def test_kirkwood_riseman_rg_over_rh_is_theta_like(protein):
    """For an ideal chain Rg/Rh from Kirkwood-Riseman sits around 1.3-1.5."""
    rg = protein.get_mean_radius_of_gyration()
    rh = protein.get_mean_hydrodynamic_radius('kirkwood-riseman')
    assert 1.2 < rg / rh < 1.6


def test_kirkwood_riseman_needs_two_residues():
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC('A').get_mean_hydrodynamic_radius('kirkwood-riseman')


def test_hydrodynamic_radius_invalid_mode(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_hydrodynamic_radius('bogus')


def test_pre_profile_shape_and_range(protein):
    idx, profile, gamma = protein.get_pre_profile(0, sample_size=200)
    n = len(protein)
    assert len(idx) == n
    assert len(profile) == n
    assert len(gamma) == n
    profile = np.asarray(profile)
    assert np.all((profile >= 0) & (profile <= np.max(profile) + 1e-9))


def test_pre_profile_boundary_rejected(protein):
    with pytest.raises(AFRCException):
        protein.get_pre_profile(len(protein))
    with pytest.raises(AFRCException):
        protein.get_pre_profile(-1)
    with pytest.raises(AFRCException):
        protein.get_pre_profile('abc')


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
def test_sampling_sizes(protein):
    assert len(protein.sample_end_to_end_distribution(n=128)) == 128
    assert len(protein.sample_radius_of_gyration_distribution(n=128)) == 128
    assert len(protein.sample_inter_residue_distance_distribution(5, 50, n=128)) == 128


def test_inter_residue_sampling_is_order_independent(protein):
    """Sampling (i, j) and (j, i) must draw from the same underlying distribution."""
    forward = protein.sample_inter_residue_distance_distribution(5, 50, n=4096)
    reverse = protein.sample_inter_residue_distance_distribution(50, 5, n=4096)
    # same distribution, so the sample means should be close
    assert np.mean(forward) == pytest.approx(np.mean(reverse), rel=0.05)


def test_inter_residue_sampling_validates_indices(protein):
    """Regression: negative indices previously wrapped round silently."""
    with pytest.raises(AFRCException):
        protein.sample_inter_residue_distance_distribution(-1, 5, n=10)
    with pytest.raises(AFRCException):
        protein.sample_inter_residue_distance_distribution(0, len(protein), n=10)


# ---------------------------------------------------------------------------
# golden-value regression snapshots
# ---------------------------------------------------------------------------
def test_golden_mean_values(protein):
    assert protein.get_mean_radius_of_gyration() == pytest.approx(30.788294820233368, abs=1e-6)
    assert protein.get_mean_radius_of_gyration('scaling law') == pytest.approx(30.80498415230296, abs=1e-6)
    assert protein.get_mean_end_to_end_distance() == pytest.approx(71.87265559462539, abs=1e-6)
    assert protein.get_mean_end_to_end_distance('distribution') == pytest.approx(71.73705605088517, abs=1e-4)
    assert protein.get_mean_hydrodynamic_radius('kirkwood-riseman') == pytest.approx(23.004186945112636, abs=1e-4)
    assert protein.get_mean_hydrodynamic_radius('nygaard') == pytest.approx(32.27795915772518, abs=1e-4)


# ---------------------------------------------------------------------------
# 0.4.3 regressions: degenerate chains, contact-fraction estimator, PRE units
# ---------------------------------------------------------------------------
def test_empty_sequence_rejected():
    """Regression: an empty sequence was accepted and gave NaN/zero results."""
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC('')


def test_hydrodynamic_radius_needs_two_residues_in_both_modes():
    """Regression: Nygaard mode returned -0.0 (with a divide-by-zero) for N = 1."""
    p = afrc.AnalyticalFRC('A')
    with pytest.raises(AFRCException):
        p.get_mean_hydrodynamic_radius('nygaard')
    with pytest.raises(AFRCException):
        p.get_mean_hydrodynamic_radius('kirkwood-riseman')


def test_contact_fraction_is_cumulative_bin_weight(protein):
    """The contact fraction is the total weight of the bins below the threshold."""
    r, p = protein.get_interresidue_distance_distribution(10, 40)
    for threshold in (5.0, 12.5, 30.0):
        expected = np.sum(p[r < threshold])
        assert protein.get_contact_fraction(10, 40, threshold) == pytest.approx(expected, rel=1e-12)
    # and a threshold beyond the grid captures everything
    assert protein.get_contact_fraction(10, 40, 1e6) == pytest.approx(1.0, rel=1e-9)


def test_contact_fraction_matches_gaussian_cdf():
    """
    Regression: the trapezoid estimator dropped half the last bin, biasing the
    fraction low by ~1.5% at a 5 A threshold. The binned cumulative sum should
    track the analytical Gaussian-chain CDF to well inside 2% even at that
    small threshold (the residual is the left-Riemann discretisation of P(r)).
    """
    from scipy.special import erf
    from afrc.polymer import PolymerObject

    p = afrc.AnalyticalFRC('A'*30)
    # the (0, 29) pair is the segment seq[0:29]; read its <r^2>
    mean_sq = PolymerObject('A'*29).RMS_Re_scaling**2
    sigma = np.sqrt(mean_sq/3.0)

    def gaussian_cdf(x):
        z = x/sigma
        return erf(z/np.sqrt(2)) - np.sqrt(2/np.pi)*z*np.exp(-z*z/2)

    for threshold, tol in ((5.0, 0.02), (10.0, 0.01), (20.0, 0.005), (40.0, 0.002)):
        assert p.get_contact_fraction(0, 29, threshold) == pytest.approx(gaussian_cdf(threshold), rel=tol)


def test_pre_profile_uses_angular_larmor_frequency(protein):
    """
    Regression: W_H defaulted to the linear proton frequency (600 MHz -> 6e8), but
    the Solomon-Bloembergen spectral-density term needs the angular frequency
    (2*pi*6e8). The default must equal the explicit angular value, and passing the
    linear value must inflate gamma by exactly the hand-computed prefactor ratio.
    """
    tau_c = 4e-9
    nu_H = 600e6

    def prefactor(omega):
        return 4*tau_c + 3*tau_c/(1 + (omega*tau_c)**2)

    np.random.seed(7)
    default_run = protein.get_pre_profile(0, sample_size=300)
    np.random.seed(7)
    angular_run = protein.get_pre_profile(0, W_H=2*np.pi*nu_H, sample_size=300)
    np.random.seed(7)
    linear_run = protein.get_pre_profile(0, W_H=nu_H, sample_size=300)

    # identical sampling, so the default and explicit-angular gamma arrays must match exactly
    assert np.array_equal(default_run[2][5], angular_run[2][5])

    # and the linear/angular ratio is the ratio of the two prefactors (~1.107 at tau_c = 4 ns)
    expected_ratio = prefactor(nu_H)/prefactor(2*np.pi*nu_H)
    assert expected_ratio > 1.05
    assert np.allclose(linear_run[2][5]/angular_run[2][5], expected_ratio, rtol=1e-9)


# ---------------------------------------------------------------------------
# inter-residue API: conventions, symmetry, calculation modes
# ---------------------------------------------------------------------------
def test_worm_like_chain_attribute_matches_sequence(all_aa):
    p = afrc.AnalyticalFRC(all_aa, adaptable_P_res=True)
    assert p.worm_like_chain.nres == len(all_aa)
    assert p.worm_like_chain.p_of_r_resolution == p.p_of_r_resolution


def test_whole_chain_accessors_delegate_to_full_sequence_polymer(protein):
    r, p = protein.get_end_to_end_distribution()
    r_po, p_po = protein.full_seq_PO.get_end_to_end_distribution()
    assert np.array_equal(r, r_po) and np.array_equal(p, p_po)


def test_interresidue_quantities_are_order_independent(protein):
    r_fwd, p_fwd = protein.get_interresidue_distance_distribution(5, 50)
    r_rev, p_rev = protein.get_interresidue_distance_distribution(50, 5)
    assert np.array_equal(r_fwd, r_rev) and np.array_equal(p_fwd, p_rev)
    assert protein.get_mean_interresidue_distance(5, 50) == protein.get_mean_interresidue_distance(50, 5)
    assert protein.get_mean_interresidue_radius_of_gyration(5, 50) == protein.get_mean_interresidue_radius_of_gyration(50, 5)
    assert protein.get_contact_fraction(5, 50, 20.0) == protein.get_contact_fraction(50, 5, 20.0)


def test_interresidue_segment_convention():
    """
    The pair (i, j) is modelled as a chain of |i - j| residues built from
    seq[i:j], so the (0, N-1) distance is a little shorter than the whole-chain
    end-to-end distance, which uses all N residues.
    """
    from afrc.polymer import PolymerObject

    seq = 'MASNDYTQQATQSYGAYPTQ'
    p = afrc.AnalyticalFRC(seq)
    n = len(seq)
    for mode in ('scaling law', 'distribution'):
        expected = PolymerObject(seq[0:n-1]).get_mean_end_to_end_distance(mode)
        assert p.get_mean_interresidue_distance(0, n-1, mode) == pytest.approx(expected, rel=1e-12)
    assert p.get_mean_interresidue_distance(0, n-1) < p.get_mean_end_to_end_distance()


def test_mean_interresidue_distance_modes_agree(protein):
    scaling = protein.get_mean_interresidue_distance(10, 60, 'scaling law')
    dist = protein.get_mean_interresidue_distance(10, 60, 'distribution')
    assert dist == pytest.approx(scaling, rel=0.005)
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_distance(10, 60, 'bogus')


def test_mean_interresidue_radius_of_gyration(protein, test_seq):
    from afrc.polymer import PolymerObject

    # a residue with itself has no extent
    assert protein.get_mean_interresidue_radius_of_gyration(7, 7) == 0.0

    # otherwise it is the Rg of the segment between the two residues
    for mode in ('scaling law', 'distribution'):
        expected = PolymerObject(test_seq[10:60]).get_mean_radius_of_gyration(mode)
        assert protein.get_mean_interresidue_radius_of_gyration(10, 60, mode) == pytest.approx(expected, rel=1e-12)

    # and the two modes agree
    scaling = protein.get_mean_interresidue_radius_of_gyration(10, 60, 'scaling law')
    dist = protein.get_mean_interresidue_radius_of_gyration(10, 60, 'distribution')
    assert dist == pytest.approx(scaling, rel=0.002)


def test_mean_interresidue_radius_of_gyration_validates_input(protein):
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_radius_of_gyration(-1, 5)
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_radius_of_gyration(0, len(protein))
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_radius_of_gyration(0, 5, 'bogus')


def test_same_residue_rg_does_not_build_the_matrix(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    assert p.get_mean_interresidue_radius_of_gyration(3, 3) == 0.0
    assert p.matrix is False


def test_inter_residue_matrix_is_built_once(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    assert p.matrix is False
    p.get_distance_map()
    first = p.matrix
    p.get_contact_map(10.0)
    assert p.matrix is first


def test_distance_map_modes_agree(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    scaling = p.get_distance_map('scaling law')
    dist = p.get_distance_map('distribution')
    upper = np.triu_indices(len(all_aa), k=1)
    assert np.allclose(dist[upper], scaling[upper], rtol=0.004)
    assert np.allclose(np.diag(dist), 0.0)


def test_internal_scaling_distribution_mode(protein):
    scaling = protein.get_internal_scaling('scaling law')
    dist = protein.get_internal_scaling('distribution')
    assert np.array_equal(scaling[:, 0], dist[:, 0])
    assert np.allclose(dist[:, 1], scaling[:, 1], rtol=0.004)


def test_same_residue_sampling_is_all_zeros(protein):
    assert np.all(protein.sample_inter_residue_distance_distribution(12, 12, n=50) == 0.0)


def test_contact_fraction_rejects_non_numeric_threshold(protein):
    """Regression: a non-numeric threshold raised a bare TypeError from NumPy."""
    with pytest.raises(AFRCException):
        protein.get_contact_fraction(10, 40, 'close')
    with pytest.raises(AFRCException):
        protein.get_contact_fraction(10, 40, None)
    # numeric strings and ints are fine
    assert protein.get_contact_fraction(10, 40, '20') == protein.get_contact_fraction(10, 40, 20)


def test_contact_map_symmetric_with_unit_diagonal(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    cm = p.get_contact_map(12.0, symmetric_map=True)
    assert np.allclose(cm, cm.T)
    assert np.allclose(np.diag(cm), 1.0)
    assert np.all((cm >= 0) & (cm <= 1))
    upper_only = p.get_contact_map(12.0)
    assert np.allclose(np.tril(upper_only, -1), 0.0)
    assert np.allclose(np.triu(upper_only), np.triu(cm))


def test_contact_fraction_decreases_with_separation(protein):
    near = protein.get_contact_fraction(20, 22, 10.0)
    far = protein.get_contact_fraction(20, 80, 10.0)
    assert near > far


def test_pre_profile_labelled_residue_and_distance_dependence(protein):
    np.random.seed(11)
    idx, profile, gamma = protein.get_pre_profile(30, sample_size=500)
    # the labelled residue is fully relaxed
    assert profile[30] == 0.0
    assert np.all(np.isinf(gamma[30]))
    # residues far from the label are less affected than those close to it
    assert profile[120] > profile[33]
    assert np.all(np.asarray(profile) <= 1.0)
    assert np.array_equal(idx, np.arange(len(protein)))


def test_nygaard_matches_published_expression(protein):
    """Nygaard et al. (2017) Eq. 7 with alpha1 = 0.216, alpha2 = 4.06, alpha3 = 0.821."""
    n = len(protein)
    rg = protein.get_mean_radius_of_gyration()
    rg_over_rh = 0.216*(rg - 4.06*n**0.33)/(n**0.60 - n**0.33) + 0.821
    assert protein.get_mean_hydrodynamic_radius('nygaard') == pytest.approx(rg/rg_over_rh, rel=1e-12)


# ---------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('bad', ['ACD EF', 'ACDEF\n', 'ACDEX', 'ACDE-F', 'ACDEFB', b'ACDEF'])
def test_malformed_sequences_raise(bad):
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC(bad)


@pytest.mark.parametrize('bad', [None, ['A', 'C'], 3.5])
def test_non_string_sequences_raise(bad):
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC(bad)


def test_mixed_case_sequence_is_accepted(all_aa):
    mixed = ''.join(c.lower() if i % 2 else c for i, c in enumerate(all_aa))
    assert afrc.AnalyticalFRC(mixed).seq == all_aa


@pytest.mark.parametrize('index', [np.int64(12), np.int32(12), '12', 12.0])
def test_residue_indices_accept_integer_like_values(protein, index):
    assert protein.get_mean_interresidue_distance(index, 40) == protein.get_mean_interresidue_distance(12, 40)


def test_last_residue_is_a_valid_index(protein):
    last = len(protein) - 1
    assert protein.get_mean_interresidue_distance(0, last) > 0
    with pytest.raises(AFRCException):
        protein.get_mean_interresidue_distance(0, last + 1)


# ---------------------------------------------------------------------------
# grid resolution
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('n', [5, 100, 400])
def test_adaptable_resolution_does_not_change_the_answer(test_seq, n):
    seq = (test_seq * 3)[:n]
    default = afrc.AnalyticalFRC(seq)
    adaptable = afrc.AnalyticalFRC(seq, adaptable_P_res=True)
    assert adaptable.get_mean_radius_of_gyration() == pytest.approx(default.get_mean_radius_of_gyration(), rel=1e-4)
    assert adaptable.get_mean_end_to_end_distance('distribution') == pytest.approx(default.get_mean_end_to_end_distance('distribution'), rel=1e-4)
    # the scaling-law values do not depend on the grid at all
    assert adaptable.get_mean_end_to_end_distance() == default.get_mean_end_to_end_distance()


def test_adaptable_resolution_is_used_for_every_distribution(all_aa):
    p = afrc.AnalyticalFRC(all_aa, adaptable_P_res=True)
    for r, _ in (p.get_end_to_end_distribution(), p.get_radius_of_gyration_distribution(),
                 p.get_interresidue_distance_distribution(2, 15)):
        assert np.allclose(np.diff(r), p.p_of_r_resolution)


# ---------------------------------------------------------------------------
# contact fraction edge cases
# ---------------------------------------------------------------------------
def test_contact_threshold_is_strict(protein):
    """A grid point exactly at the threshold is not counted."""
    r, p = protein.get_interresidue_distance_distribution(10, 30)
    threshold = r[200]
    assert protein.get_contact_fraction(10, 30, threshold) == pytest.approx(np.sum(p[:200]), rel=1e-12)


def test_contact_fraction_limits(protein):
    assert protein.get_contact_fraction(10, 30, -5.0) == 0.0
    assert protein.get_contact_fraction(10, 30, 0.0) == 0.0
    assert protein.get_contact_fraction(10, 30, 1e9) == pytest.approx(1.0, rel=1e-12)


def test_contact_map_rejects_non_numeric_threshold(all_aa):
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC(all_aa).get_contact_map('far')


def test_distance_map_distribution_mode_symmetric(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    dm = p.get_distance_map('distribution', symmetric_map=True)
    assert np.allclose(dm, dm.T)
    assert np.allclose(np.triu(dm), p.get_distance_map('distribution'))


# ---------------------------------------------------------------------------
# sampling statistics
# ---------------------------------------------------------------------------
def test_sampling_is_reproducible_with_a_seed(protein):
    np.random.seed(17)
    a = (protein.sample_end_to_end_distribution(100), protein.sample_radius_of_gyration_distribution(100),
         protein.sample_inter_residue_distance_distribution(3, 60, 100))
    np.random.seed(17)
    b = (protein.sample_end_to_end_distribution(100), protein.sample_radius_of_gyration_distribution(100),
         protein.sample_inter_residue_distance_distribution(3, 60, 100))
    for x, y in zip(a, b):
        assert np.array_equal(x, y)


def test_sample_means_match_the_distributions(protein):
    np.random.seed(23)
    n = 40000
    re = protein.sample_end_to_end_distribution(n)
    rg = protein.sample_radius_of_gyration_distribution(n)
    pair = protein.sample_inter_residue_distance_distribution(3, 60, n)
    assert np.mean(re) == pytest.approx(protein.get_mean_end_to_end_distance('distribution'), rel=0.01)
    assert np.mean(rg) == pytest.approx(protein.get_mean_radius_of_gyration('distribution'), rel=0.01)
    assert np.mean(pair) == pytest.approx(protein.get_mean_interresidue_distance(3, 60, 'distribution'), rel=0.01)


def test_samples_lie_on_the_distribution_grid(protein):
    r, _ = protein.get_interresidue_distance_distribution(3, 60)
    assert np.all(np.isin(protein.sample_inter_residue_distance_distribution(3, 60, 300), r))


# ---------------------------------------------------------------------------
# PRE profile parameters
# ---------------------------------------------------------------------------
def _pre(protein, **kwargs):
    np.random.seed(31)
    return protein.get_pre_profile(40, sample_size=400, **kwargs)


def test_pre_profile_is_reproducible_with_a_seed(protein):
    a, b = _pre(protein), _pre(protein)
    assert np.array_equal(a[1], b[1])
    assert all(np.array_equal(x, y) for x, y in zip(a[2], b[2]))


def test_pre_gamma_has_one_value_per_sample(protein):
    _, _, gamma = _pre(protein)
    assert all(len(g) == 400 for g in gamma)


def test_pre_gamma_scales_with_the_spectral_density(protein):
    """Changing tau_c rescales every Gamma_2 by the ratio of the spectral-density prefactors."""
    omega = 2 * np.pi * 600e6

    def prefactor(tau_c_ns):
        tau = tau_c_ns * 1e-9
        return 4 * tau + 3 * tau / (1 + (omega * tau)**2)

    slow, fast = _pre(protein, tau_c=8), _pre(protein, tau_c=2)
    expected = prefactor(8) / prefactor(2)
    for i in (0, 20, 100):
        assert np.allclose(slow[2][i] / fast[2][i], expected, rtol=1e-12)


def test_pre_profile_responds_to_experimental_parameters(protein):
    """Longer delays and correlation times give more relaxation; a larger R_2D gives less."""
    base = np.asarray(_pre(protein)[1])
    others = np.arange(len(base)) != 40
    assert np.all(np.asarray(_pre(protein, t_delay=30)[1])[others] <= base[others])
    assert np.all(np.asarray(_pre(protein, tau_c=10)[1])[others] <= base[others])
    assert np.all(np.asarray(_pre(protein, R_2D=40)[1])[others] >= base[others])


def test_pre_profile_is_near_one_far_from_the_label(protein):
    profile = np.asarray(_pre(protein)[1])
    assert profile[153] > 0.8
    assert profile[41] < 0.1
