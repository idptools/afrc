"""
Tests for ``afrc/polymer_models/nudep_saw.py`` - the nu-dependent SAW model.
"""

import numpy as np
import pytest

from afrc.polymer_models import nudep_saw


def test_distribution_is_valid_pmf(all_aa):
    model = nudep_saw.NuDepSAW(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_and_rms_distinct_and_ordered(all_aa):
    model = nudep_saw.NuDepSAW(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert mean > 0
    assert rms >= mean


def test_mean_radius_of_gyration_positive(all_aa):
    assert nudep_saw.NuDepSAW(all_aa).get_mean_radius_of_gyration() > 0


def test_larger_nu_expands_chain(all_aa):
    """A larger scaling exponent should give a larger end-to-end distance."""
    model = nudep_saw.NuDepSAW(all_aa)
    compact = model.get_mean_end_to_end_distance(nu=0.4)
    expanded = model.get_mean_end_to_end_distance(nu=0.6)
    assert expanded > compact


def test_sampling_size(all_aa):
    model = nudep_saw.NuDepSAW(all_aa)
    assert len(model.sample_end_to_end_distribution(n=100)) == 100


def test_exception_is_real_exception():
    assert issubclass(nudep_saw.NuDepSAWException, Exception)
