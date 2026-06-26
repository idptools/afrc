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


def test_exception_is_real_exception():
    assert issubclass(wlc2.WLC2Exception, Exception)


def test_rejects_short_sequence():
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2('AA')


def test_rejects_non_positive_lp(all_aa):
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2(all_aa, lp=-1)
