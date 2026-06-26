"""
Tests for ``afrc/polymer_models/frc.py`` - the freely rotating chain model.
"""

import numpy as np
import pytest

from afrc.polymer_models import frc


def test_distribution_is_valid_pmf(all_aa):
    model = frc.FreelyRotatingChain(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_and_rms_distinct_and_ordered(all_aa):
    model = frc.FreelyRotatingChain(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert mean > 0
    assert rms >= mean


def test_rms_matches_analytic_characteristic_ratio(all_aa):
    """RMS should match sqrt(<R^2>) from the finite-N freely-rotating-chain formula."""
    b, c_inf = 3.8, 2.0
    N = len(all_aa)
    model = frc.FreelyRotatingChain(all_aa, b=b, c_inf=c_inf)
    rms = model.get_root_mean_squared_end_to_end_distance()

    alpha = (c_inf - 1.0) / (c_inf + 1.0)
    mean_sq = c_inf * N * b**2 - 2 * b**2 * alpha * (1 - alpha**N) / (1 - alpha)**2
    assert rms == pytest.approx(np.sqrt(mean_sq), rel=1e-3)


def test_c_inf_one_recovers_ideal_chain(all_aa):
    """With c_inf = 1 (freely jointed limit) the RMS is exactly b*sqrt(N)."""
    b = 3.8
    model = frc.FreelyRotatingChain(all_aa, b=b, c_inf=1.0)
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rms == pytest.approx(b * np.sqrt(len(all_aa)), rel=1e-3)


def test_larger_c_inf_expands_chain(all_aa):
    flexible = frc.FreelyRotatingChain(all_aa, c_inf=1.0).get_root_mean_squared_end_to_end_distance()
    stiff = frc.FreelyRotatingChain(all_aa, c_inf=4.0).get_root_mean_squared_end_to_end_distance()
    assert stiff > flexible


def test_longer_bond_expands_chain(all_aa):
    small = frc.FreelyRotatingChain(all_aa, b=3.0).get_mean_end_to_end_distance()
    large = frc.FreelyRotatingChain(all_aa, b=4.5).get_mean_end_to_end_distance()
    assert large > small


def test_mean_radius_of_gyration_relationship(all_aa):
    model = frc.FreelyRotatingChain(all_aa)
    rg = model.get_mean_radius_of_gyration()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rg == pytest.approx(rms / np.sqrt(6))


def test_sampling_size(all_aa):
    model = frc.FreelyRotatingChain(all_aa)
    assert len(model.sample_end_to_end_distribution(n=100)) == 100


def test_exception_is_real_exception():
    assert issubclass(frc.FRCException, Exception)


def test_rejects_non_positive_bond_length(all_aa):
    with pytest.raises(frc.FRCException):
        frc.FreelyRotatingChain(all_aa, b=0)


def test_rejects_non_positive_c_inf(all_aa):
    with pytest.raises(frc.FRCException):
        frc.FreelyRotatingChain(all_aa, c_inf=0)
