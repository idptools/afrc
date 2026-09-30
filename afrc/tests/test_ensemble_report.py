"""
Tests for ``afrc/ensemble_report.py`` - checking a generated ensemble against its model.

All ensembles are drawn with fixed seeds, so these tests are deterministic.
"""

import numpy as np
import pytest

import afrc
from afrc.afrc import AFRCException
from afrc.ensemble import freely_rotating_chain_msd, mean_squared_distance_map, sample_freely_rotating_chain
from afrc.ensemble_report import (KS_P_THRESHOLD, MIN_CONFORMATIONS, Z_THRESHOLD, Check, EnsembleReport, ModelExpectations,
                                  _separation_bands, _z_scores, compare_ensemble_to_gaussian_model, compare_ensemble_to_model)
from afrc.tests.conftest import TEST_SEQ

P53 = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'


def _report(seq, n, seed, transform=None, **kwargs):
    p = afrc.AnalyticalFRC(seq)
    xyz = p.sample_conformations(n, seed=seed)
    if transform is not None:
        xyz = transform(xyz)
    return p, compare_ensemble_to_gaussian_model(xyz, p.get_mean_squared_distance_map(), **kwargs)


# ---------------------------------------------------------------------------
# a correct ensemble passes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('seq, n', [(P53, 1000), (TEST_SEQ, 300), ('GS', 500), ('ACDEF', 200), ('W' * 40, 2000)])
def test_correct_ensembles_pass(seq, n):
    _, report = _report(seq, n, seed=11)
    assert report.assessed
    assert report.passed
    assert all(check.passed for check in report.checks)
    assert 'Verdict: PASS' in report.format()


def test_no_false_alarms_across_many_seeds():
    """A calibration guard: correct ensembles should essentially never be flagged."""
    p = afrc.AnalyticalFRC(P53[:20])
    D = p.get_mean_squared_distance_map()
    failures = sum(not compare_ensemble_to_gaussian_model(p.sample_conformations(200, seed=s), D).passed for s in range(60))
    assert failures == 0


# ---------------------------------------------------------------------------
# a wrong ensemble is flagged
# ---------------------------------------------------------------------------
def _failed(report):
    return {check.name.split(' (')[0] for check in report.checks if check.passed is False}


@pytest.mark.parametrize('scale, n', [(1.02, 1000), (0.98, 1000), (1.01, 5000)])
@pytest.mark.parametrize('seed', range(3))
def test_small_size_errors_are_flagged(scale, n, seed):
    """2% errors are always caught with 1000 conformations, and 1% errors with 5000."""
    _, report = _report(P53, n, seed=seed, transform=lambda xyz: scale * xyz)
    assert not report.passed
    assert 'Mean inter-residue distances' in _failed(report)
    assert 'Verdict: WARNING' in report.format()


def test_three_percent_size_error_fails_several_checks():
    _, report = _report(P53, 1000, seed=1, transform=lambda xyz: 1.03 * xyz)
    assert {'Kirkwood-Riseman hydrodynamic radius', 'Mean-squared inter-residue distances',
            'Mean inter-residue distances'} <= _failed(report)


def test_ensemble_from_the_wrong_sequence_is_flagged():
    p = afrc.AnalyticalFRC(P53)
    wrong = afrc.AnalyticalFRC('G' * len(P53)).sample_conformations(1000, seed=1)
    assert not compare_ensemble_to_gaussian_model(wrong, p.get_mean_squared_distance_map()).passed


def test_shuffled_beads_are_flagged():
    rng = np.random.default_rng(0)
    _, report = _report(P53, 1000, seed=1, transform=lambda xyz: xyz[:, rng.permutation(len(P53))])
    assert {'Mean first-to-last bead distance', 'Inter-residue distance distributions'} <= _failed(report)


def test_fixed_bond_random_walk_is_flagged():
    """A freely jointed chain scaled to the right overall size still has the wrong local structure."""
    rng = np.random.default_rng(0)
    steps = rng.standard_normal((1000, len(P53) - 1, 3))
    steps *= 3.8 / np.linalg.norm(steps, axis=-1, keepdims=True)
    walk = np.concatenate([np.zeros((1000, 1, 3)), np.cumsum(steps, axis=1)], axis=1)
    p = afrc.AnalyticalFRC(P53)
    scale = np.sqrt(np.mean(p.get_mean_squared_distance_map()[0, 1:] / (3.8**2 * np.arange(1, len(P53)))))
    report = compare_ensemble_to_gaussian_model(scale * walk, p.get_mean_squared_distance_map())
    assert 'Inter-residue distance distributions' in _failed(report)
    assert not report.passed


# ---------------------------------------------------------------------------
# the model values are the AFRC's own
# ---------------------------------------------------------------------------
def test_model_values_match_the_afrc_methods():
    p, report = _report(P53, 200, seed=3)
    by_name = {check.name.split(' (')[0]: check.values for check in report.checks}
    n = len(P53)
    assert by_name['Kirkwood-Riseman hydrodynamic radius']['model'] == pytest.approx(p.get_mean_hydrodynamic_radius(), rel=1e-10)
    assert by_name['Mean first-to-last bead distance']['model'] == pytest.approx(
        p.get_mean_interresidue_distance(0, n - 1, 'distribution'), rel=1e-4)
    # <Rg^2> = (1/N^2) sum_{i<j} <r_ij^2>
    D = p.get_mean_squared_distance_map()
    assert by_name['Root-mean-square radius of gyration']['model'] == pytest.approx(np.sqrt(np.sum(np.triu(D)) / n**2), rel=1e-12)


def test_context_rows():
    p, report = _report(P53, 500, seed=3, reference_mean_rg=17.7, reference_mean_re=41.3)
    names = [check.name for check in report.context]
    assert any('<Rg>' in name for name in names)
    assert any('whole-chain <Re>' in name for name in names)
    assert any('Adjacent bead spacing' in name for name in names)
    assert any('scaling exponent' in name for name in names)
    assert all(check.passed is None for check in report.context)
    nu = next(check for check in report.context if 'scaling exponent' in check.name).values
    assert nu['ensemble'] == pytest.approx(nu['model'], abs=0.02)


def test_context_rows_are_optional():
    _, report = _report(P53, 200, seed=3)
    assert not any('<Rg>' in check.name for check in report.context)


# ---------------------------------------------------------------------------
# edge cases and structure
# ---------------------------------------------------------------------------
def test_single_bead_is_not_assessed():
    _, report = _report('A', 500, seed=0)
    assert report.checks == []
    assert not report.assessed
    assert not report.passed
    assert 'NOT ASSESSED' in report.format()


def test_too_few_conformations_are_not_assessed():
    _, report = _report(P53, MIN_CONFORMATIONS - 1, seed=0)
    assert not report.assessed
    assert all(check.passed is None for check in report.checks)
    assert 'NOT ASSESSED' in report.format()
    # ...but the numbers are still reported
    assert 'ensemble' in report.format()


def test_two_residues_are_assessed():
    _, report = _report('GS', 500, seed=0)
    assert report.passed
    assert not any('scaling exponent' in check.name for check in report.context)


def test_pair_budget_limits_the_pair_based_checks():
    _, report = _report(P53, 500, seed=0, max_pair_evaluations=1225 * 150)
    assert report.n_conformations == 500
    assert report.n_pair_conformations == 150
    assert 'pair-based checks use the first 150' in report.format()
    assert report.passed


def test_mismatched_inputs_raise():
    p = afrc.AnalyticalFRC(P53)
    D = p.get_mean_squared_distance_map()
    with pytest.raises(AFRCException):
        compare_ensemble_to_gaussian_model(np.zeros((10, 5, 3)), D)
    with pytest.raises(AFRCException):
        compare_ensemble_to_gaussian_model(np.zeros((10, len(P53))), D)
    with pytest.raises(AFRCException):
        compare_ensemble_to_gaussian_model(np.zeros((0, len(P53), 3)), D)


@pytest.mark.parametrize('n_res, expected', [
    (2, [(1, 1)]),
    (4, [(1, 1), (2, 3)]),
    (50, [(1, 1), (2, 3), (4, 10), (11, 30), (31, 49)]),
    (2000, [(1, 1), (2, 3), (4, 10), (11, 30), (31, 100), (101, 300), (301, 1000), (1001, 1999)]),
])
def test_separation_bands_cover_every_separation_once(n_res, expected):
    bands = _separation_bands(n_res)
    assert bands == expected
    covered = [k for low, high in bands for k in range(low, high + 1)]
    assert covered == list(range(1, n_res))


def test_report_formatting():
    report = EnsembleReport('Toy', 3, 100, 100,
                            checks=[Check('A check', True, 'line one\nline two', {'z': 0.1}),
                                    Check('Another', False, 'bad', {'z': 9.0})],
                            context=[Check('Context', None, 'fyi')])
    text = report.format(['Header line'])
    assert text.splitlines()[:3] == ['Toy ensemble report', '===================', 'Header line']
    assert '[OK]    A check' in text and '        line two' in text
    assert '[FAIL]  Another' in text and '[info]  Context' in text
    assert 'Verdict: WARNING - 1 check(s) outside the expected range: Another.' in text
    assert f'({Z_THRESHOLD:g} standard errors)' in text


def test_thresholds_are_sensible():
    assert Z_THRESHOLD >= 4
    assert KS_P_THRESHOLD <= 1e-4
    assert MIN_CONFORMATIONS >= 50


# ---------------------------------------------------------------------------
# general models: structural checks, context and notes
# ---------------------------------------------------------------------------
def _frc(n_res=20, n=300, cos_angle=1/3, seed=0):
    xyz = sample_freely_rotating_chain(n_res, 3.8, cos_angle, n, np.random.default_rng(seed))
    D = mean_squared_distance_map(n_res, lambda k: freely_rotating_chain_msd(k, 3.8, cos_angle))
    return xyz, D


def test_structural_checks_pass_for_the_right_geometry():
    xyz, D = _frc()
    report = compare_ensemble_to_model(xyz, ModelExpectations('FRC', D, bond_length=3.8, bond_angle_cosine=1/3,
                                                              contour_length_per_residue=3.8))
    structural = {check.name: check for check in report.checks if check.structural}
    assert set(structural) == {'Bond length', 'Bond angle', 'Finite extensibility'}
    assert all(check.passed for check in structural.values())
    assert report.passed
    assert '109.47 degrees' in structural['Bond angle'].summary


@pytest.mark.parametrize('field, value', [('bond_length', 3.9), ('bond_angle_cosine', 0.3), ('contour_length_per_residue', 3.0)])
def test_structural_checks_fail_for_the_wrong_geometry(field, value):
    xyz, D = _frc()
    kwargs = {'bond_length': 3.8, 'bond_angle_cosine': 1/3, 'contour_length_per_residue': 3.8, field: value}
    report = compare_ensemble_to_model(xyz, ModelExpectations('FRC', D, **kwargs))
    assert not report.passed
    assert 'Verdict: WARNING' in report.format()


def test_structural_failures_are_reported_even_with_few_conformations():
    xyz, D = _frc(n=10)
    report = compare_ensemble_to_model(xyz, ModelExpectations('FRC', D, bond_length=4.0))
    assert not report.assessed
    assert report.failed == ['Bond length']
    assert 'Verdict: WARNING' in report.format()


def test_non_gaussian_models_skip_the_gaussian_checks():
    xyz, D = _frc()
    names = [check.name for check in compare_ensemble_to_model(xyz, ModelExpectations('FRC', D)).checks]
    assert 'Root-mean-square first-to-last bead distance' in names
    assert not any('Kirkwood-Riseman' in name or 'Kolmogorov' in name or 'Mean inter-residue' in name for name in names)


def test_end_to_end_distribution_context_and_notes():
    xyz, D = _frc(n=500)
    r = np.arange(0, 100, 0.05)
    p = r**2 * np.exp(-3 * r**2 / (2 * D[0, -1]))
    report = compare_ensemble_to_model(xyz, ModelExpectations('FRC', D, end_to_end_distribution=(r, p / p.sum()),
                                                              end_to_end_label='a Gaussian', notes=('a caveat',)))
    row = next(check for check in report.context if 'distribution vs a Gaussian' in check.name)
    assert 0 <= row.values['cdf_difference'] < 0.1
    assert report.notes == ['a caveat']
    assert '* a caveat' in report.format()


def test_z_scores_for_quantities_that_do_not_fluctuate():
    expected = np.array([1.0, 1.0, 2.0])
    # rounding-level agreement with a rounding-level standard error: no alarm
    z = _z_scores(np.array([1.0 + 2e-16, 1.0, 2.0]), expected, np.array([1e-17, 0.0, 1e-3]))
    assert np.all(np.abs(z) < 1)
    # a real mismatch in a fixed quantity is flagged loudly
    assert abs(_z_scores(np.array([1.000001]), np.array([1.0]), np.array([0.0]))[0]) > 100
    # an unknown standard error stays unknown
    assert np.isnan(_z_scores(np.array([1.0]), np.array([1.0]), np.array([np.nan]))[0])
