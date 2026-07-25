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


def test_rg_uses_the_same_nu_as_the_size_scaling(all_aa):
    """Regression: the Rg ratio hard-coded nu = 0.589 while Ree scaled as N^0.598.

    Both now read self.nu, so changing that single attribute must move Rg and Re
    together.
    """
    model = saw.SAW(all_aa)
    N = len(all_aa)

    gamma, nu = model.gamma, model.nu
    ratio = gamma*(gamma + 1) / (2*(gamma + 2*nu)*(gamma + 2*nu + 1))
    expected = model.get_root_mean_squared_end_to_end_distance() * np.sqrt(ratio)
    assert model.get_mean_radius_of_gyration() == pytest.approx(expected, rel=1e-12)

    # and nu really is the exponent driving the size scaling
    assert model.get_root_mean_squared_end_to_end_distance() == pytest.approx(5.5 * N**nu, rel=1e-3)


def test_rg_matches_nu_dependent_saw(all_aa):
    """With a consistent nu, SAW and NuDepSAW(nu=0.598) give the same Rg."""
    from afrc.polymer_models import nudep_saw

    fixed = saw.SAW(all_aa)
    variable = nudep_saw.NuDepSAW(all_aa)
    assert fixed.get_mean_radius_of_gyration() == pytest.approx(
        variable.get_mean_radius_of_gyration(nu=fixed.nu), rel=1e-4)


def test_prefactor_scales_dimensions(all_aa):
    model = saw.SAW(all_aa)
    small = model.get_mean_end_to_end_distance(prefactor=4.0)
    large = model.get_mean_end_to_end_distance(prefactor=6.0)
    assert large > small


def test_rms_equals_prefactor_times_n_to_the_nu(all_aa):
    """sqrt(<Re^2>) should come out as prefactor * N^0.598 by construction."""
    N = len(all_aa)
    model = saw.SAW(all_aa)
    for prefactor in (4.0, 5.5, 6.5):
        rms = model.get_root_mean_squared_end_to_end_distance(prefactor=prefactor)
        assert rms == pytest.approx(prefactor * N**0.598, rel=1e-3)


def test_long_chain_grid_covers_the_distribution():
    """Regression: for long chains the r-grid used to truncate the tail of P(r).

    Ree scales as N^0.598 while the default grid scaled as N^0.5, so the grid
    eventually cut into the distribution and biased the mean/RMS low.
    """
    N = 2000
    model = saw.SAW('A'*N)
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rms == pytest.approx(5.5 * N**0.598, rel=1e-3)

    # and essentially none of the probability mass should sit at the grid edge
    dist, prob = model.get_end_to_end_distribution()
    assert prob[-1] < 1e-12


def test_exception_is_real_exception():
    assert issubclass(saw.SAWException, Exception)
