"""
Tests for ``afrc/polymer_models/wlc.py`` - the Zhou worm-like chain model.
"""

import numpy as np
import pytest

from afrc.polymer_models import wlc


def test_distribution_is_valid_pmf(all_aa):
    model = wlc.WormLikeChain(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_no_negative_probabilities(all_aa):
    """Regression: the zeta series can go negative in the tail; values are clamped."""
    model = wlc.WormLikeChain(all_aa)
    _, prob = model.get_end_to_end_distribution()
    assert np.all(prob >= 0)


def test_mean_and_rms_positive_and_ordered(all_aa):
    model = wlc.WormLikeChain(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert mean > 0
    assert rms >= mean


def test_exception_is_real_exception():
    assert issubclass(wlc.WLCException, Exception)


def test_rejects_non_positive_lp(all_aa):
    with pytest.raises(wlc.WLCException):
        wlc.WormLikeChain(all_aa, lp=0)


def test_rejects_non_positive_aa_size(all_aa):
    with pytest.raises(wlc.WLCException):
        wlc.WormLikeChain(all_aa, aa_size=0)


def test_no_probability_beyond_contour_length():
    """Regression: the Zhou series assigned real weight to r > Lc for short chains.

    A chain can never be longer than its contour length, but for a 2-residue chain
    at lp = 3 A roughly 30% of the distribution (and over half at lp = 10 A for a
    5-residue chain) previously sat beyond Lc.
    """
    b = 3.8
    for n, lp in ((2, 3.0), (5, 3.0), (5, 6.0), (10, 10.0), (50, 3.0)):
        model = wlc.WormLikeChain('A' * n, lp=lp, aa_size=b)
        dist, prob = model.get_end_to_end_distribution()
        assert np.all(np.isfinite(prob))
        assert np.sum(prob[dist > n * b]) == 0.0
        assert np.sum(prob) == pytest.approx(1.0)


def test_rms_matches_exact_worm_like_chain():
    """For Lc >> Lp the distribution reproduces the exact WLC second moment.

    <R^2> = 2 Lp Lc - 2 Lp^2 (1 - exp(-Lc/Lp)); this also pins the sign convention
    of the zeta correction series.
    """
    b = 3.8
    for n, lp in ((50, 3.0), (100, 6.0), (300, 10.0)):
        Lc = n * b
        exact = np.sqrt(2 * lp * Lc - 2 * lp**2 * (1 - np.exp(-Lc / lp)))
        rms = wlc.WormLikeChain('A' * n, lp=lp, aa_size=b).get_root_mean_squared_end_to_end_distance()
        assert rms == pytest.approx(exact, rel=1e-4)


def test_grid_covers_tail_for_stiff_chains():
    """Regression: the fixed 21*sqrt(N) grid truncated the tail of P(r) for large lp."""
    dist, prob = wlc.WormLikeChain('A' * 100, lp=10.0).get_end_to_end_distribution()
    assert prob[-1] < 1e-8


def test_chain_shorter_than_persistence_length_raises():
    """With Lc < lp the expansion has no valid region; this must not return NaNs."""
    model = wlc.WormLikeChain('AA', lp=10.0)
    with pytest.raises(wlc.WLCException):
        model.get_end_to_end_distribution()
