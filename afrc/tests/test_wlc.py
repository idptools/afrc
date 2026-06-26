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
