"""
Interface-contract tests shared by every polymer model in the package.

Each model is built by a small factory so the same checks run against all seven:
valid distributions on a regular grid, means and RMS values that are consistent
with those distributions, sensible dependence on chain length and grid spacing,
reproducible sampling, and pinned regression values.
"""

import numpy as np
import pytest

import afrc
from afrc.polymer_models.fjc import FreelyJointedChain
from afrc.polymer_models.frc import FreelyRotatingChain
from afrc.polymer_models.nudep_saw import NuDepSAW
from afrc.polymer_models.saw import SAW
from afrc.polymer_models.wlc import WormLikeChain
from afrc.polymer_models.wlc2 import WormLikeChain2


class Model:
    """Uniform wrapper so the contract tests read the same for every model."""

    def __init__(self, name, build, mean_re, has_rg=True, can_sample=True, custom_resolution=True):
        self.name = name
        self.build = build
        self.mean_re = mean_re
        self.has_rg = has_rg
        self.can_sample = can_sample
        self.custom_resolution = custom_resolution

    def __repr__(self):
        return self.name


def _afrc(seq, res=None):
    return afrc.AnalyticalFRC(seq)


MODELS = [
    # the AFRC's default mean is the scaling law; the distribution mode is the
    # one that must agree with its own P(r)
    Model('AnalyticalFRC', _afrc, lambda m: m.get_mean_end_to_end_distance('distribution'), custom_resolution=False),
    Model('FreelyJointedChain', lambda s, res=0.05: FreelyJointedChain(s, res), lambda m: m.get_mean_end_to_end_distance()),
    Model('FreelyRotatingChain', lambda s, res=0.05: FreelyRotatingChain(s, res), lambda m: m.get_mean_end_to_end_distance()),
    Model('WormLikeChain', lambda s, res=0.05: WormLikeChain(s, res), lambda m: m.get_mean_end_to_end_distance(), has_rg=False, can_sample=False),
    Model('WormLikeChain2', lambda s, res=0.05: WormLikeChain2(s, res), lambda m: m.get_mean_end_to_end_distance(), can_sample=False),
    Model('SAW', lambda s, res=0.05: SAW(s, res), lambda m: m.get_mean_end_to_end_distance(), can_sample=False),
    Model('NuDepSAW', lambda s, res=0.05: NuDepSAW(s, res), lambda m: m.get_mean_end_to_end_distance()),
]

SEQ = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQ'


def _rms(model):
    r, p = model.get_end_to_end_distribution()
    return np.sqrt(np.sum(p * r**2))


# ---------------------------------------------------------------------------
# distributions
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('model', MODELS, ids=repr)
def test_distribution_is_a_valid_pmf_on_a_regular_grid(model):
    r, p = model.build(SEQ).get_end_to_end_distribution()
    assert r.shape == p.shape
    assert r.dtype.kind == 'f' and p.dtype.kind == 'f'
    assert r[0] == 0.0
    assert np.allclose(np.diff(r), 0.05)
    assert np.all(np.isfinite(p))
    assert np.all(p >= 0)
    assert np.sum(p) == pytest.approx(1.0, abs=1e-12)
    # nothing sits at exactly zero distance for a non-degenerate chain
    assert p[0] == 0.0


@pytest.mark.parametrize('model', [m for m in MODELS if m.custom_resolution], ids=repr)
def test_custom_resolution_is_honoured_and_does_not_change_the_answer(model):
    fine = model.build(SEQ, 0.05)
    coarse = model.build(SEQ, 0.5)
    r, _ = coarse.get_end_to_end_distribution()
    assert np.allclose(np.diff(r), 0.5)
    assert model.mean_re(coarse) == pytest.approx(model.mean_re(fine), rel=0.005)


@pytest.mark.parametrize('model', MODELS, ids=repr)
def test_distribution_is_deterministic(model):
    a = model.build(SEQ).get_end_to_end_distribution()
    b = model.build(SEQ).get_end_to_end_distribution()
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


# ---------------------------------------------------------------------------
# means and root-mean-square values
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('model', MODELS, ids=repr)
def test_mean_is_the_expectation_of_the_distribution(model):
    m = model.build(SEQ)
    r, p = m.get_end_to_end_distribution()
    assert model.mean_re(m) == pytest.approx(np.sum(r * p), rel=1e-12)


@pytest.mark.parametrize('model', [m for m in MODELS if m.name != 'AnalyticalFRC'], ids=repr)
def test_rms_is_consistent_with_the_distribution(model):
    m = model.build(SEQ)
    assert m.get_root_mean_squared_end_to_end_distance() == pytest.approx(_rms(m), rel=1e-12)


@pytest.mark.parametrize('model', MODELS, ids=repr)
def test_rms_exceeds_mean(model):
    """Jensen's inequality: sqrt(<r^2>) >= <r>, strictly for a spread-out distribution."""
    m = model.build(SEQ)
    assert _rms(m) > model.mean_re(m) > 0


@pytest.mark.parametrize('model', MODELS, ids=repr)
def test_longer_chains_are_larger(model):
    means = [model.mean_re(model.build('A' * n)) for n in (10, 20, 40, 80, 160)]
    assert np.all(np.diff(means) > 0)


@pytest.mark.parametrize('model', [m for m in MODELS if m.has_rg], ids=repr)
def test_radius_of_gyration_is_smaller_than_the_end_to_end_distance(model):
    m = model.build(SEQ)
    rg = m.get_mean_radius_of_gyration()
    assert 0 < rg < model.mean_re(m)
    # and for an ideal-ish chain it is roughly Re/sqrt(6); the excluded-volume
    # models sit a little lower, the stiff ones a little higher
    assert 0.3 < rg / _rms(m) < 0.5


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('model', [m for m in MODELS if m.can_sample], ids=repr)
def test_sampling_is_reproducible_and_on_the_grid(model):
    m = model.build(SEQ)
    np.random.seed(21)
    a = m.sample_end_to_end_distribution(n=500)
    np.random.seed(21)
    b = m.sample_end_to_end_distribution(n=500)
    assert np.array_equal(a, b)
    assert len(a) == 500

    r, _ = m.get_end_to_end_distribution()
    assert np.all(np.isin(a, r))


@pytest.mark.parametrize('model', [m for m in MODELS if m.can_sample], ids=repr)
def test_sample_moments_match_the_distribution(model):
    m = model.build(SEQ)
    r, p = m.get_end_to_end_distribution()
    mean = np.sum(r * p)
    sd = np.sqrt(np.sum(p * (r - mean)**2))

    np.random.seed(8)
    samples = m.sample_end_to_end_distribution(n=40000)
    # the standard error of the mean is sd/200, so 5 standard errors is very safe
    assert np.mean(samples) == pytest.approx(mean, abs=5 * sd / 200)
    assert np.std(samples) == pytest.approx(sd, rel=0.03)


# ---------------------------------------------------------------------------
# pinned regression values (defaults, 60-residue sequence)
# ---------------------------------------------------------------------------
GOLDEN = {
    'WormLikeChain': (33.93358794852771, 36.74234613130608, None),
    'WormLikeChain2': (34.06112046310616, 36.81187484011443, 14.80654334278507),
    'FreelyJointedChain': (26.91523725586497, 29.189402719882544, 11.916523760053282),
    'FreelyRotatingChain': (38.111162943892644, 41.365927939599075, 16.88756936479267),
    'SAW': (59.7232227012581, 63.635139439947395, 25.341677974602845),
    'NuDepSAW': (39.55621848281252, 42.60281679875327, 18.259392774434396),
}


@pytest.mark.parametrize('model', [m for m in MODELS if m.name in GOLDEN], ids=repr)
def test_golden_values(model):
    mean, rms, rg = GOLDEN[model.name]
    m = model.build(SEQ)
    assert m.get_mean_end_to_end_distance() == pytest.approx(mean, rel=1e-9)
    assert m.get_root_mean_squared_end_to_end_distance() == pytest.approx(rms, rel=1e-9)
    if rg is not None:
        assert m.get_mean_radius_of_gyration() == pytest.approx(rg, rel=1e-9)
