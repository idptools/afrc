"""
Tests for ``afrc/polymer_models/fjc.py`` - the freely jointed chain model.
"""

import numpy as np
import pytest

from afrc.polymer_models import fjc


def test_distribution_is_valid_pmf(all_aa):
    model = fjc.FreelyJointedChain(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_finite_extensibility(all_aa):
    """No probability mass should sit beyond the contour length L = N*b."""
    model = fjc.FreelyJointedChain(all_aa, b=3.8)
    dist, prob = model.get_end_to_end_distribution()
    L = len(all_aa) * 3.8
    assert np.all(dist < L)
    assert np.all(prob[dist >= L] == 0) if np.any(dist >= L) else True


def test_mean_and_rms_distinct_and_ordered(all_aa):
    model = fjc.FreelyJointedChain(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert mean > 0
    assert rms >= mean


def test_rms_close_to_ideal_chain(all_aa):
    """In the long-chain limit the RMS should approach b*sqrt(N)."""
    b = 3.8
    model = fjc.FreelyJointedChain(all_aa, b=b)
    rms = model.get_root_mean_squared_end_to_end_distance()
    ideal = b * np.sqrt(len(all_aa))
    # finite extensibility pulls the RMS slightly below the ideal value
    assert rms == pytest.approx(ideal, rel=0.1)
    assert rms <= ideal * 1.001


def test_mean_radius_of_gyration_relationship(all_aa):
    model = fjc.FreelyJointedChain(all_aa)
    rg = model.get_mean_radius_of_gyration()
    rms = model.get_root_mean_squared_end_to_end_distance()
    assert rg == pytest.approx(rms / np.sqrt(6))


def test_longer_segment_expands_chain(all_aa):
    small = fjc.FreelyJointedChain(all_aa, b=3.0).get_mean_end_to_end_distance()
    large = fjc.FreelyJointedChain(all_aa, b=4.5).get_mean_end_to_end_distance()
    assert large > small


def test_sampling_size(all_aa):
    model = fjc.FreelyJointedChain(all_aa)
    assert len(model.sample_end_to_end_distribution(n=100)) == 100


def test_exception_is_real_exception():
    assert issubclass(fjc.FJCException, Exception)


def test_rejects_non_positive_segment_length(all_aa):
    with pytest.raises(fjc.FJCException):
        fjc.FreelyJointedChain(all_aa, b=0)


def test_zero_length_chain_is_degenerate():
    model = fjc.FreelyJointedChain('')
    dist, prob = model.get_end_to_end_distribution()
    assert dist == pytest.approx([0.0])
    assert prob == pytest.approx([1.0])
    assert np.all(model.sample_end_to_end_distribution(n=20) == 0.0)


@pytest.mark.parametrize('b', [3.8, 10.0, 20.0])
def test_grid_covers_tail_for_large_segments(b):
    """
    Regression: the grid was fixed at 21*sqrt(N), which truncated the tail for
    segment lengths above ~7 A (the RMS came out over 25% low at b = 20 A).
    """
    n = 50
    model = fjc.FreelyJointedChain('A' * n, b=b)
    dist, prob = model.get_end_to_end_distribution()
    assert prob[-1] < 1e-8 * prob.max()
    rms = model.get_root_mean_squared_end_to_end_distance()
    # finite extensibility pulls the RMS slightly below the ideal b*sqrt(N)
    assert rms < b * np.sqrt(n)
    assert rms == pytest.approx(b * np.sqrt(n), rel=0.02)


def test_grid_never_extends_past_contour_length():
    model = fjc.FreelyJointedChain('A' * 5, b=20.0)
    dist, _ = model.get_end_to_end_distribution()
    assert dist[-1] < 5 * 20.0


def test_sampling_reuses_the_cached_distribution(all_aa):
    np.random.seed(2)
    model = fjc.FreelyJointedChain(all_aa * 3)
    mean = model.get_mean_end_to_end_distance()
    samples = model.sample_end_to_end_distribution(n=4000)
    assert np.mean(samples) == pytest.approx(mean, rel=0.05)
