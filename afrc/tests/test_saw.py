"""
Tests for ``afrc/polymer_models/saw.py`` - the self-avoiding walk model.
"""

import numpy as np
import pytest

from afrc.polymer_models import saw


def test_distribution_is_valid_pmf(all_aa):
    model = saw.SAW(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_and_rms_distinct_and_ordered(all_aa):
    model = saw.SAW(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    # RMS is always >= mean for a non-degenerate distribution
    assert mean > 0
    assert rms >= mean


def test_mean_radius_of_gyration_positive(all_aa):
    assert saw.SAW(all_aa).get_mean_radius_of_gyration() > 0


def test_prefactor_scales_dimensions(all_aa):
    model = saw.SAW(all_aa)
    small = model.get_mean_end_to_end_distance(prefactor=4.0)
    large = model.get_mean_end_to_end_distance(prefactor=6.0)
    assert large > small


def test_exception_is_real_exception():
    assert issubclass(saw.SAWException, Exception)
