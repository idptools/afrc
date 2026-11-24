"""
Tests for 3D ensemble generation from the auxiliary polymer models: the freely
jointed chain, freely rotating chain, both worm-like chains, the SAW and the
nu-dependent SAW.

A shared contract runs against every model (shape, seeding, validation,
writing, and a passing report for a correct ensemble), followed by
model-specific tests of the exact quantities each generator must reproduce and
of the report's ability to catch a wrong ensemble.
"""

import os

import numpy as np
import pytest

from afrc.afrc import AFRCException
from afrc.ensemble import (discrete_worm_like_chain_msd, freely_rotating_chain_msd, sample_worm_like_chain,
                           worm_like_chain_discretization, worm_like_chain_msd)
from afrc.polymer_models.fjc import FreelyJointedChain
from afrc.polymer_models.frc import FreelyRotatingChain
from afrc.polymer_models.nudep_saw import NuDepSAW, NuDepSAWException
from afrc.polymer_models.saw import SAW, SAWException
from afrc.polymer_models.wlc import WormLikeChain
from afrc.polymer_models.wlc2 import WLC2Exception, WormLikeChain2

P53 = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'


class Case:
    """A model plus the per-call parameters its ensemble methods take."""

    def __init__(self, name, build, **kwargs):
        self.name = name
        self.build = build
        self.kwargs = kwargs

    def __repr__(self):
        return self.name


CASES = [
    Case('fjc', lambda s: FreelyJointedChain(s)),
    Case('fjc b=5', lambda s: FreelyJointedChain(s, b=5.0)),
    Case('frc', lambda s: FreelyRotatingChain(s)),
    Case('frc c_inf=0.5', lambda s: FreelyRotatingChain(s, c_inf=0.5)),
    Case('frc c_inf=9', lambda s: FreelyRotatingChain(s, c_inf=9.0)),
    Case('wlc', lambda s: WormLikeChain(s)),
    Case('wlc lp=20', lambda s: WormLikeChain(s, lp=20.0)),
    Case('wlc lp=0.5', lambda s: WormLikeChain(s, lp=0.5)),
    Case('wlc2', lambda s: WormLikeChain2(s)),
    Case('saw', lambda s: SAW(s)),
    Case('saw prefactor=6', lambda s: SAW(s), prefactor=6.0),
    Case('saw-nu nu=0.35', lambda s: NuDepSAW(s), nu=0.35),
    Case('saw-nu nu=0.588', lambda s: NuDepSAW(s), nu=0.588, prefactor=6.0),
]


def _sample(case, seq=P53, n=200, seed=1):
    model = case.build(seq)
    return model, model.sample_conformations(n, seed=seed, **case.kwargs)


# ---------------------------------------------------------------------------
# the shared contract
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('case', CASES, ids=repr)
def test_shape_and_centring(case):
    _, xyz = _sample(case, n=40)
    assert xyz.shape == (40, len(P53), 3)
    assert np.all(np.isfinite(xyz))
    assert np.allclose(xyz.mean(axis=1), 0.0, atol=1e-9)


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_seeding(case):
    model = case.build(P53)
    a = model.sample_conformations(10, seed=5, **case.kwargs)
    assert np.array_equal(a, model.sample_conformations(10, seed=5, **case.kwargs))
    assert np.array_equal(a, model.sample_conformations(10, seed=np.random.default_rng(5), **case.kwargs))
    assert not np.array_equal(a, model.sample_conformations(10, seed=6, **case.kwargs))
    np.random.seed(3)
    expected = np.random.random()
    np.random.seed(3)
    model.sample_conformations(10, seed=5, **case.kwargs)
    assert np.random.random() == expected


@pytest.mark.parametrize('case', CASES, ids=repr)
@pytest.mark.parametrize('bad', [0, -1, 2.5, 'ten', None])
def test_number_of_conformations_is_validated(case, bad):
    with pytest.raises(AFRCException):
        case.build(P53).sample_conformations(bad, **case.kwargs)


@pytest.mark.parametrize('case', CASES, ids=repr)
@pytest.mark.parametrize('seq', ['', 'A'])
def test_degenerate_chains(case, seq):
    if case.name == 'wlc2' and seq == '':
        pytest.skip("the O'Brien model rejects chains shorter than lp, including empty ones")
    xyz = case.build(seq).sample_conformations(4, seed=0, **case.kwargs)
    assert xyz.shape == (4, len(seq), 3)
    assert np.all(xyz == 0.0)


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_correct_ensembles_pass(case):
    model, xyz = _sample(case, n=1000, seed=11)
    report = model.check_ensemble(xyz, **case.kwargs)
    assert report.passed, report.format()
    assert 'Verdict: PASS' in report.format()


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_first_to_last_mean_square_is_the_model_value(case):
    model, xyz = _sample(case, n=4000, seed=12)
    msd = model.get_mean_squared_distance_map(**case.kwargs)
    first_to_last = np.mean(np.sum((xyz[:, -1] - xyz[:, 0])**2, axis=-1))
    # the relative standard error of <r^2> is ~1.3% here
    assert first_to_last == pytest.approx(msd[0, -1], rel=0.07)


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_mean_squared_distance_map_structure(case):
    msd = case.build(P53).get_mean_squared_distance_map(**case.kwargs)
    assert msd.shape == (len(P53), len(P53))
    assert np.array_equal(msd, msd.T)
    assert np.all(np.diag(msd) == 0)
    # homogeneous chains: depends only on |i - j|, and grows with it
    assert np.allclose(np.diagonal(msd, offset=7), msd[0, 7])
    assert np.all(np.diff(msd[0]) > 0)


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_save_ensemble_pair(case, tmp_path):
    md = pytest.importorskip('mdtraj')
    model = case.build(P53.lower())
    xyz = model.save_ensemble(str(tmp_path / 'ens'), n=15, seed=2, **case.kwargs)
    assert sorted(os.listdir(tmp_path)) == ['ens.pdb', 'ens.xtc']
    traj = md.load(str(tmp_path / 'ens.xtc'), top=str(tmp_path / 'ens.pdb'))
    assert traj.n_frames == 15
    assert ''.join(residue.code for residue in traj.topology.residues) == P53
    assert np.allclose(traj.xyz * 10, xyz, atol=0.006)
    assert np.array_equal(xyz, model.sample_conformations(15, seed=2, **case.kwargs))


@pytest.mark.parametrize('case', CASES, ids=repr)
def test_save_ensemble_pdb_only(case, tmp_path):
    model = case.build(P53)
    model.save_ensemble(str(tmp_path / 'ens.pdb'), n=4, seed=2, pdb_only=True, **case.kwargs)
    assert os.listdir(tmp_path) == ['ens.pdb']
    with open(tmp_path / 'ens.pdb') as fh:
        text = fh.read()
    assert text.count('MODEL ') == 4
    assert 'ensemble' in text.splitlines()[0]


def test_saving_needs_a_real_sequence(tmp_path):
    with pytest.raises(AFRCException, match='non-standard'):
        FreelyJointedChain('ACDXZ').save_ensemble(str(tmp_path / 'bad'), n=2, pdb_only=True)
    assert os.listdir(tmp_path) == []


# ---------------------------------------------------------------------------
# freely jointed chain
# ---------------------------------------------------------------------------
def test_fjc_mean_squared_distances_are_k_b_squared():
    msd = FreelyJointedChain(P53, b=4.2).get_mean_squared_distance_map()
    k = np.abs(np.subtract.outer(np.arange(len(P53)), np.arange(len(P53))))
    assert np.allclose(msd, 4.2**2 * k, rtol=1e-12)


def test_fjc_bonds_are_exact_and_directions_are_uncorrelated():
    _, xyz = _sample(Case('fjc', lambda s: FreelyJointedChain(s)), n=5000, seed=3)
    bonds = np.diff(xyz, axis=1)
    assert np.allclose(np.linalg.norm(bonds, axis=-1), 3.8, rtol=1e-12)
    cosines = np.sum(bonds[:, 1:] * bonds[:, :-1], axis=-1) / 3.8**2
    # uniformly random directions: <cos> = 0 and <cos^2> = 1/3
    assert np.mean(cosines) == pytest.approx(0.0, abs=0.01)
    assert np.mean(cosines**2) == pytest.approx(1 / 3, abs=0.01)


# ---------------------------------------------------------------------------
# freely rotating chain
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('c_inf', [0.5, 2.0, 9.0])
def test_frc_mean_squared_distances_match_bond_correlations(c_inf):
    """<r^2> = sum over bond pairs of b^2 alpha^|p - q|, computed directly."""
    b, n = 3.8, 12
    alpha = (c_inf - 1) / (c_inf + 1)
    msd = FreelyRotatingChain('A' * n, b=b, c_inf=c_inf).get_mean_squared_distance_map()
    for k in range(1, n):
        p = np.arange(k)
        direct = b * b * np.sum(alpha ** np.abs(np.subtract.outer(p, p)))
        assert msd[0, k] == pytest.approx(direct, rel=1e-12)


def test_frc_whole_chain_formula_matches_the_model():
    """The pair map uses the same exact expression as the model's own P(r)."""
    model = FreelyRotatingChain('A' * 31, c_inf=3.0)
    msd = model.get_mean_squared_distance_map()
    r, p = FreelyRotatingChain('A' * 30, c_inf=3.0, p_of_r_resolution=0.01).get_end_to_end_distribution()
    assert msd[0, -1] == pytest.approx(np.sum(r * r * p), rel=1e-4)


@pytest.mark.parametrize('c_inf, angle', [(2.0, 109.4712), (1.0, 90.0), (9.0, 143.1301)])
def test_frc_bond_angles_are_exact(c_inf, angle):
    _, xyz = _sample(Case('frc', lambda s: FreelyRotatingChain(s, c_inf=c_inf)), n=100, seed=4)
    bonds = np.diff(xyz, axis=1)
    units = bonds / np.linalg.norm(bonds, axis=-1, keepdims=True)
    cosines = np.sum(units[:, 1:] * units[:, :-1], axis=-1)
    assert np.allclose(180 - np.degrees(np.arccos(cosines)), angle, atol=1e-4)
    assert np.allclose(np.linalg.norm(bonds, axis=-1), 3.8, rtol=1e-12)


def test_frc_torsions_are_uniform():
    _, xyz = _sample(Case('frc', lambda s: FreelyRotatingChain(s)), n=3000, seed=5)
    b1, b2, b3 = (xyz[:, i + 1] - xyz[:, i] for i in (10, 11, 12))
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    torsion = np.arctan2(np.sum(np.cross(n1, n2) * b2, axis=-1) / np.linalg.norm(b2, axis=-1), np.sum(n1 * n2, axis=-1))
    counts, _ = np.histogram(torsion, bins=12, range=(-np.pi, np.pi))
    assert np.all(np.abs(counts - 250) < 5 * np.sqrt(250))


# ---------------------------------------------------------------------------
# worm-like chains
# ---------------------------------------------------------------------------
def test_wlc_mean_squared_distances_are_the_exact_worm_like_chain():
    model = WormLikeChain(P53, lp=4.0, aa_size=3.6)
    k = np.arange(1, len(P53))
    assert np.allclose(model.get_mean_squared_distance_map()[0, 1:], worm_like_chain_msd(3.6 * k, 4.0), rtol=1e-12)
    assert np.array_equal(model.get_mean_squared_distance_map(), WormLikeChain2(P53, lp=4.0, aa_size=3.6).get_mean_squared_distance_map())


def test_wlc_ensembles_do_not_depend_on_the_analytical_form():
    """Zhou and O'Brien describe the same chain, so they generate the same ensemble."""
    assert np.array_equal(WormLikeChain(P53).sample_conformations(20, seed=9), WormLikeChain2(P53).sample_conformations(20, seed=9))


def test_wlc_discretization_is_cached_and_the_note_reports_it():
    model = WormLikeChain(P53)
    assert model._discretization is None
    m, correlation = model._ensemble_discretization()
    assert model._discretization == (m, correlation)
    assert m > 1 and 0 < correlation < 1
    report = model.check_ensemble(model.sample_conformations(150, seed=1))
    assert any(f'{m} straight sub-segments' in note for note in report.notes)


# ---------------------------------------------------------------------------
# worm-like chain discretization
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('lp', [0.1, 0.5, 1.0, 3.0, 10.0, 50.0])
def test_discretization_meets_the_tolerance(lp):
    n_res = 134
    m, correlation = worm_like_chain_discretization(3.8, lp, n_res)
    k = np.arange(1, n_res, dtype=float)
    deviation = discrete_worm_like_chain_msd(k, 3.8, m, correlation) / worm_like_chain_msd(3.8 * k, lp) - 1
    assert np.max(np.abs(deviation)) <= 2e-4
    # and one fewer sub-segment would not (i.e. m is the smallest that works for this rule)
    # the two rules can agree to ~1e-5, so compare exactly rather than with isclose
    exponential = abs(correlation - np.exp(-3.8 / m / lp)) < 1e-14
    fewer = np.exp(-3.8 / (m - 1) / lp) if exponential else (2 * lp - 3.8 / (m - 1)) / (2 * lp + 3.8 / (m - 1))
    if m > 1 and 0 < fewer < 1:
        assert np.max(np.abs(discrete_worm_like_chain_msd(k, 3.8, m - 1, fewer) / worm_like_chain_msd(3.8 * k, lp) - 1)) > 2e-4


def test_flexible_chains_use_the_matched_correlation():
    """Regression: lp = 0.1 A needed 807 sub-segments per residue (and ~12 s for 500 conformations)."""
    m, correlation = worm_like_chain_discretization(3.8, 0.1, 134)
    s = 3.8 / m
    assert correlation == pytest.approx((2 * 0.1 - s) / (2 * 0.1 + s), rel=1e-12)
    assert m <= 250


def test_stiff_chains_keep_the_worm_like_chain_correlation():
    m, correlation = worm_like_chain_discretization(3.8, 50.0, 134)
    assert correlation == pytest.approx(np.exp(-3.8 / m / 50.0), rel=1e-12)


def test_matched_correlation_reproduces_the_long_range_size():
    """s(1 + c)/(1 - c) = 2 lp, so <r^2> grows at exactly the continuous chain's rate."""
    m, c = worm_like_chain_discretization(3.8, 0.2, 60)
    s = 3.8 / m
    assert s * (1 + c) / (1 - c) == pytest.approx(2 * 0.2, rel=1e-12)


def test_discrete_msd_is_a_freely_rotating_chain_of_sub_segments():
    k = np.arange(1, 20, dtype=float)
    assert np.allclose(discrete_worm_like_chain_msd(k, 3.8, 7, 0.6), freely_rotating_chain_msd(7 * k, 3.8 / 7, 0.6), rtol=1e-14)


def test_sampler_reproduces_the_flexible_discretization_exactly():
    n_res = 12
    m, c = worm_like_chain_discretization(3.8, 0.1, n_res)
    xyz = sample_worm_like_chain(n_res, 3.8, 20000, np.random.default_rng(4), m, c)
    for k in (1, 4, 11):
        observed = np.mean(np.sum((xyz[:, k:] - xyz[:, :-k])**2, axis=-1))
        assert observed == pytest.approx(discrete_worm_like_chain_msd(np.array([k]), 3.8, m, c)[0], rel=0.02)


def test_sub_segment_directions_follow_the_von_mises_fisher_distribution():
    """With one sub-segment per residue the bonds are the sub-segments themselves."""
    from afrc.ensemble import _mean_cosine_to_concentration

    xyz = sample_worm_like_chain(200, 2.0, 1000, np.random.default_rng(5), 1, 0.7)
    bonds = np.diff(xyz, axis=1)
    assert np.allclose(np.linalg.norm(bonds, axis=-1), 2.0, rtol=1e-12)
    cosines = np.sum(bonds[:, 1:] * bonds[:, :-1], axis=-1) / 4.0
    kappa = _mean_cosine_to_concentration(0.7)
    assert np.mean(cosines) == pytest.approx(0.7, abs=0.003)
    # for a von Mises-Fisher distribution in 3D, <cos^2> = 1 - 2 <cos> / kappa
    assert np.mean(cosines**2) == pytest.approx(1 - 2 * 0.7 / kappa, abs=0.003)


def test_block_wise_random_draws_are_exact():
    """Very large ensembles draw random numbers in blocks of sub-segments; the chain must not change."""
    n_conf = 300_000          # makes the block (1M / n_conf = 3) smaller than the 10 sub-segments
    xyz = sample_worm_like_chain(3, 3.8, n_conf, np.random.default_rng(6), 10, 0.8)
    observed = np.mean(np.sum((xyz[:, 2] - xyz[:, 0])**2, axis=-1))
    assert observed == pytest.approx(discrete_worm_like_chain_msd(np.array([2]), 3.8, 10, 0.8)[0], rel=0.005)


@pytest.mark.parametrize('bad', [0.0, 1.0, -0.2, 1.5])
def test_sampler_rejects_invalid_correlations(bad):
    with pytest.raises(AFRCException):
        sample_worm_like_chain(5, 3.8, 10, np.random.default_rng(0), 4, bad)


def test_wlc_is_finitely_extensible():
    model, xyz = _sample(Case('wlc lp=50', lambda s: WormLikeChain(s, lp=50.0)), n=300, seed=6)
    k = np.abs(np.subtract.outer(np.arange(len(P53)), np.arange(len(P53))))
    distances = np.linalg.norm(xyz[:, :, None] - xyz[:, None], axis=-1)
    mask = k > 0
    assert np.all(distances[:, mask] <= 3.8 * k[mask] * (1 + 1e-9))
    # a stiff chain gets close to its contour length
    assert np.max(distances[:, 0, -1]) > 0.8 * 3.8 * (len(P53) - 1)


def test_wlc_reports_for_chains_too_short_for_the_analytical_forms():
    """The context rows are skipped when the Zhou or O'Brien forms cannot be evaluated."""
    zhou = WormLikeChain('AA', lp=10.0)
    report = zhou.check_ensemble(zhou.sample_conformations(200, seed=1))
    assert report.passed
    assert not any('distance distribution' in check.name for check in report.context)

    obrien = WormLikeChain2('AA', lp=5.0)    # N - 1 = 1 residue is shorter than lp
    report = obrien.check_ensemble(obrien.sample_conformations(200, seed=1))
    assert report.passed
    assert not any('distance distribution' in check.name for check in report.context)


def test_wlc2_still_rejects_chains_shorter_than_lp():
    with pytest.raises(WLC2Exception):
        WormLikeChain2('A', lp=5.0)


# ---------------------------------------------------------------------------
# self-avoiding walks (Gaussian approximation)
# ---------------------------------------------------------------------------
def test_saw_mean_squared_distances_follow_the_scaling_law():
    msd = SAW(P53).get_mean_squared_distance_map(prefactor=6.0)
    k = np.arange(1, len(P53))
    assert np.allclose(msd[0, 1:], (6.0 * k**0.598)**2, rtol=1e-12)
    msd = NuDepSAW(P53).get_mean_squared_distance_map(nu=0.4, prefactor=5.0)
    assert np.allclose(msd[0, 1:], (5.0 * k**0.4)**2, rtol=1e-12)


def test_saw_ensembles_scale_with_the_prefactor():
    model = SAW(P53)
    base = model.sample_conformations(10, seed=3, prefactor=1.0)
    assert np.allclose(model.sample_conformations(10, seed=3, prefactor=6.5), 6.5 * base, rtol=1e-12)
    # the covariance factor is computed once and reused
    assert len(model._unit_factors) == 1


def test_saw_nu_caches_one_factor_per_exponent():
    model = NuDepSAW(P53)
    for nu in (0.4, 0.5, 0.4, 0.5):
        model.sample_conformations(5, seed=0, nu=nu)
    assert set(model._unit_factors) == {0.4, 0.5}


def test_saw_nu_at_0598_matches_the_fixed_exponent_saw():
    assert np.allclose(NuDepSAW(P53).sample_conformations(10, seed=4, nu=0.598, prefactor=5.5),
                       SAW(P53).sample_conformations(10, seed=4, prefactor=5.5), rtol=1e-10)


def test_saw_reports_are_marked_as_a_gaussian_approximation():
    model, xyz = _sample(Case('saw', lambda s: SAW(s)), n=500, seed=7)
    report = model.check_ensemble(xyz)
    assert any('Gaussian approximation' in note for note in report.notes)
    assert any("des Cloizeaux" in check.name for check in report.context)
    assert 'Gaussian approximation' in report.format()


@pytest.mark.parametrize('call', [
    lambda: SAW(P53).sample_conformations(5, prefactor=0),
    lambda: SAW(P53).get_mean_squared_distance_map(prefactor=-1),
    lambda: SAW(P53).check_ensemble(np.zeros((5, len(P53), 3)), prefactor=0),
])
def test_saw_parameters_are_validated(call):
    with pytest.raises(SAWException):
        call()


@pytest.mark.parametrize('kwargs', [{'nu': 0}, {'nu': 1}, {'nu': 1.5}, {'prefactor': 0}, {'prefactor': -2}])
def test_saw_nu_parameters_are_validated(kwargs):
    with pytest.raises(NuDepSAWException):
        NuDepSAW(P53).sample_conformations(5, **kwargs)
    with pytest.raises(NuDepSAWException):
        NuDepSAW(P53).get_mean_squared_distance_map(**kwargs)


# ---------------------------------------------------------------------------
# the report catches wrong ensembles
# ---------------------------------------------------------------------------
def _failed(report):
    return {check.name.split(' (')[0] for check in report.checks if check.passed is False}


@pytest.mark.parametrize('model, wrong, kwargs, expected', [
    (FreelyJointedChain(P53), lambda: FreelyJointedChain(P53, b=3.81).sample_conformations(1000, seed=1), {},
     {'Bond length', 'Mean-squared inter-residue distances'}),
    (FreelyRotatingChain(P53), lambda: FreelyRotatingChain(P53, c_inf=2.05).sample_conformations(1000, seed=1), {},
     {'Bond angle', 'Mean-squared inter-residue distances'}),
    (FreelyRotatingChain(P53, c_inf=1.0), lambda: FreelyJointedChain(P53).sample_conformations(1000, seed=1), {},
     {'Bond angle'}),
    (WormLikeChain(P53), lambda: WormLikeChain(P53, lp=3.1).sample_conformations(1000, seed=1), {},
     {'Mean-squared inter-residue distances'}),
    (SAW(P53), lambda: SAW(P53).sample_conformations(1000, seed=1, prefactor=5.61), {},
     {'Mean-squared inter-residue distances'}),
    (SAW(P53), lambda: NuDepSAW(P53).sample_conformations(1000, seed=1, nu=0.5), {},
     {'Root-mean-square radius of gyration', 'Mean-squared inter-residue distances'}),
])
def test_wrong_ensembles_are_flagged(model, wrong, kwargs, expected):
    report = model.check_ensemble(wrong(), **kwargs)
    assert not report.passed
    assert expected <= _failed(report)


def test_a_stretched_chain_breaks_finite_extensibility():
    model = WormLikeChain(P53, lp=50.0)
    xyz = 1.3 * model.sample_conformations(300, seed=2)
    assert 'Finite extensibility' in _failed(model.check_ensemble(xyz))


@pytest.mark.parametrize('case', [c for c in CASES if c.name != 'wlc lp=20'], ids=repr)
def test_single_residue_report(case):
    """One bead has nothing to check; the report says so rather than failing."""
    model = case.build('W')
    report = model.check_ensemble(model.sample_conformations(150, seed=0, **case.kwargs), **case.kwargs)
    assert report.checks == []
    assert not report.assessed
    assert 'NOT ASSESSED' in report.format()


@pytest.mark.parametrize('kwargs', [{'nu': 1.2}, {'prefactor': -1}])
def test_saw_nu_check_ensemble_validates_parameters(kwargs):
    with pytest.raises(NuDepSAWException):
        NuDepSAW(P53).check_ensemble(np.zeros((5, len(P53), 3)), **kwargs)


@pytest.mark.parametrize('bad', [None, 42, ''])
def test_validate_sequence_rejects_non_sequences(bad):
    from afrc.ensemble import validate_sequence

    with pytest.raises(AFRCException):
        validate_sequence(bad)


def test_validate_sequence_uppercases():
    from afrc.ensemble import validate_sequence

    assert validate_sequence('acdEF') == 'ACDEF'
