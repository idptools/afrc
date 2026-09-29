"""
Tests for ``afrc/polymer_models/wlc2.py`` - the O'Brien worm-like chain model.
"""

import numpy as np
import pytest

from afrc.polymer_models import wlc2


def test_distribution_is_valid_pmf(all_aa):
    model = wlc2.WormLikeChain2(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_rms_and_rg_positive(all_aa):
    model = wlc2.WormLikeChain2(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    rg = model.get_mean_radius_of_gyration()
    assert mean > 0
    assert rms >= mean
    assert rg > 0


def test_rg_matches_benoit_doty(all_aa):
    """Regression: the <Rg^2> expression previously had the wrong sign on the -lp^2 term."""
    lp, b = 3.0, 3.8
    N = len(all_aa)
    Lc = N * b

    model = wlc2.WormLikeChain2(all_aa, lp=lp, aa_size=b)

    # standard Benoit-Doty worm-like chain result
    rg_sq = Lc*lp/3 - lp**2 + 2*lp**3/Lc - 2*lp**4/Lc**2*(1 - np.exp(-Lc/lp))
    assert model.get_mean_radius_of_gyration() == pytest.approx(np.sqrt(rg_sq), rel=1e-9)


def test_rg_never_exceeds_the_rigid_rod_bound():
    """A chain can be no more extended than a rigid rod, so Rg^2 <= Lc^2/12.

    Regression: with the sign error on the -lp^2 term this bound was violated for
    short chains (e.g. a 2-residue chain gave Rg^2 = 21.1 against a rod bound of 4.8).
    """
    b = 3.8
    for n in (2, 3, 5, 10, 50, 200):
        seq = 'A' * n
        Lc = n * b
        rg = wlc2.WormLikeChain2(seq, lp=3.0, aa_size=b).get_mean_radius_of_gyration()
        assert rg**2 <= Lc**2/12


def test_rg_grows_with_persistence_length(all_aa):
    flexible = wlc2.WormLikeChain2(all_aa, lp=2.0).get_mean_radius_of_gyration()
    stiff = wlc2.WormLikeChain2(all_aa, lp=5.0).get_mean_radius_of_gyration()
    assert stiff > flexible


def test_exception_is_real_exception():
    assert issubclass(wlc2.WLC2Exception, Exception)


def test_rejects_chain_shorter_than_persistence_length():
    """The contour length (N * aa_size), not the residue count, is compared with lp."""
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2('AA', lp=20.0)

    # 2 residues is a contour length of 7.6 A, comfortably above the default lp = 3 A
    assert wlc2.WormLikeChain2('AA').get_mean_radius_of_gyration() > 0


def test_rejects_non_positive_lp(all_aa):
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2(all_aa, lp=-1)


def test_no_probability_at_or_beyond_contour_length():
    """The O'Brien expression is undefined for r >= Lc; those points carry no weight."""
    n, b = 5, 3.8
    dist, prob = wlc2.WormLikeChain2('A' * n, aa_size=b).get_end_to_end_distribution()
    assert np.all(np.isfinite(prob))
    assert np.sum(prob[dist >= n * b]) == 0.0
    assert np.sum(prob) == pytest.approx(1.0)


def test_rejects_non_positive_aa_size(all_aa):
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2(all_aa, aa_size=0)


def test_empty_sequence_rejected():
    """An empty sequence has zero contour length, so fails the contour-length check."""
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2('')


def _exact_wlc_rms(lp, lc):
    return np.sqrt(2*lp*lc - 2*lp**2*(1 - np.exp(-lc/lp)))


@pytest.mark.parametrize('n, lp', [(600, 2.0), (800, 3.0), (1000, 3.0), (5000, 3.0)])
def test_long_chains_do_not_overflow(n, lp):
    """
    Regression: C1 ~ exp(3*Lc/(4*Lp)) overflowed once that exponent passed ~709
    (about 500 residues at lp = 2 A, 750 at lp = 3 A) and the distribution came
    back all NaN.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter('error')
        model = wlc2.WormLikeChain2('A' * n, lp=lp)
        dist, prob = model.get_end_to_end_distribution()

    assert np.all(np.isfinite(prob))
    assert np.sum(prob) == pytest.approx(1.0)
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rms == pytest.approx(_exact_wlc_rms(lp, n * 3.8), rel=0.002)


@pytest.mark.parametrize('lp, n, aa_size', [(50.0, 1000, 3.8), (3.0, 100, 10.0), (30.0, 10, 3.8)])
def test_grid_covers_tail_for_stiff_chains_and_large_segments(lp, n, aa_size):
    """
    Regression: the grid depended on lp but not aa_size and was capped, so stiff
    chains and large segment sizes had their tail cut off (RMS ~2-3% low).
    """
    model = wlc2.WormLikeChain2('A' * n, lp=lp, aa_size=aa_size)
    dist, prob = model.get_end_to_end_distribution()
    assert prob[-1] < 1e-8 * prob.max()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rms == pytest.approx(_exact_wlc_rms(lp, n * aa_size), rel=0.01)


def test_grid_never_extends_past_contour_length():
    model = wlc2.WormLikeChain2('A' * 10, lp=30.0)
    dist, _ = model.get_end_to_end_distribution()
    assert dist[-1] < 10 * 3.8


def test_c1_matches_closed_form_normalization():
    """C1 is now built in log space; it must still equal O'Brien's expression."""
    model = wlc2.WormLikeChain2('A' * 20)
    a = model.alpha
    expected = 1.0/(np.pi**1.5 * np.exp(-a) * a**-1.5 * (1 + 3/a + 15/(4*a**2)))
    assert model.C1 == pytest.approx(expected, rel=1e-12)


def test_analytic_normalization_matches_numerical():
    """With C1 and the Lc^3 denominator the O'Brien P(r) integrates to one, and the log-space form is the same function."""
    model = wlc2.WormLikeChain2('A' * 40, p_of_r_resolution=0.01)
    lc = 40 * 3.8
    dist, prob = model.get_end_to_end_distribution()
    x2 = (dist / lc)**2
    analytic = 4*np.pi*model.C1*dist**2 / (lc**3*(1 - x2)**4.5) * np.exp(-model.alpha/(1 - x2))
    assert np.sum(analytic) * 0.01 == pytest.approx(1.0, rel=1e-3)
    assert np.allclose(prob, analytic/np.sum(analytic), rtol=1e-9, atol=1e-15)
