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


def test_rms_equals_prefactor_times_n_to_the_nu(all_aa):
    """The model is normalised so that sqrt(<Re^2>) is exactly prefactor * N^nu.

    Regression: the A1/A2 gamma-function prefactors were previously computed with the
    wrong arguments/exponents, which broke this identity and was being papered over
    with an unexplained factor of pi in the size scale.
    """
    N = len(all_aa)
    model = nudep_saw.NuDepSAW(all_aa)
    for nu in (0.35, 0.45, 0.5, 0.55, 0.598):
        for prefactor in (4.0, 5.5):
            rms = model.get_root_mean_squared_end_to_end_distance(nu=nu, prefactor=prefactor)
            assert rms == pytest.approx(prefactor * N**nu, rel=1e-3)


def test_matches_fixed_exponent_saw(all_aa):
    """At nu = 0.598 this should reproduce the fixed-exponent SAW model."""
    from afrc.polymer_models import saw

    fixed = saw.SAW(all_aa).get_mean_end_to_end_distance()
    variable = nudep_saw.NuDepSAW(all_aa).get_mean_end_to_end_distance(nu=0.598)
    assert variable == pytest.approx(fixed, rel=0.01)


def test_rejects_out_of_range_nu(all_aa):
    """nu must sit strictly between 0 and 1 (nu = 1 previously divided by zero)."""
    model = nudep_saw.NuDepSAW(all_aa)
    for bad_nu in (0.0, 1.0, -0.2, 1.5):
        with pytest.raises(nudep_saw.NuDepSAWException):
            model.get_mean_end_to_end_distance(nu=bad_nu)


def test_exception_is_real_exception():
    assert issubclass(nudep_saw.NuDepSAWException, Exception)


def test_zero_length_chain_is_degenerate():
    model = nudep_saw.NuDepSAW('')
    dist, prob = model.get_end_to_end_distribution(nu=0.5)
    assert dist == pytest.approx([0.0])
    assert prob == pytest.approx([1.0])


@pytest.mark.parametrize('prefactor', [0, -1.0])
def test_rejects_non_positive_prefactor(all_aa, prefactor):
    """Regression: a non-positive prefactor silently returned NaN."""
    model = nudep_saw.NuDepSAW(all_aa)
    with pytest.raises(nudep_saw.NuDepSAWException):
        model.get_mean_end_to_end_distance(prefactor=prefactor)
    with pytest.raises(nudep_saw.NuDepSAWException):
        model.sample_end_to_end_distribution(n=10, prefactor=prefactor)


@pytest.mark.parametrize('nu', [0.35, 0.5, 0.588])
def test_a1_normalizes_the_distribution(nu):
    """A1 should make the analytical P(r) integrate to one for any nu."""
    from scipy.integrate import quad

    model = nudep_saw.NuDepSAW('A')
    g = (model.gamma - 1)/nu
    delta = 1/(1 - nu)
    a1 = model._NuDepSAW__compute_A1(delta, g)
    a2 = model._NuDepSAW__compute_A2(delta, g)
    integral, _ = quad(lambda x: 4*np.pi*a1 * x**(2 + g) * np.exp(-a2 * x**delta), 0, np.inf)
    assert integral == pytest.approx(1.0, rel=1e-8)


def test_sampling_respects_nu(all_aa):
    np.random.seed(3)
    model = nudep_saw.NuDepSAW(all_aa * 5)
    compact = model.sample_end_to_end_distribution(n=4000, nu=0.4)
    expanded = model.sample_end_to_end_distribution(n=4000, nu=0.6)
    assert np.mean(expanded) > np.mean(compact)
