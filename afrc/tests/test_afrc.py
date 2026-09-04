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

    p = afrc.AnalyticalFRC('A'*30)
    # force the inter-residue matrix to exist, then read the pair's <r^2>
    p.get_interresidue_distance_distribution(0, 29)
    mean_sq = p.matrix[0][29].RMS_Re_scaling**2
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
