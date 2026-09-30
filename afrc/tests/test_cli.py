"""
Tests for the ``afrc-ensemble`` command-line tool (``afrc/cli.py``).

The command is exercised through ``ensemble_main(argv)``, which is exactly what the
installed console script calls.
"""

import os
import re
import subprocess
import sys

import numpy as np
import pytest

import afrc
from afrc.cli import ENSEMBLE_MODELS, build_ensemble_parser, ensemble_main

P53 = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'


@pytest.fixture(autouse=True)
def _run_in_a_temporary_directory(tmp_path, monkeypatch):
    """
    Run every test from its own temporary directory, so a command that uses the
    default output prefix can never write into the repository.
    """
    monkeypatch.chdir(tmp_path)


def _run(tmp_path, *args):
    """Run afrc-ensemble writing into tmp_path; returns the exit status."""
    return ensemble_main(['-s', P53, '-o', str(tmp_path / 'out'), *args])


def _read(path):
    with open(path) as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------------------
def test_defaults():
    args = build_ensemble_parser().parse_args(['-s', 'ACDE'])
    assert args.number_of_conformers == 1000
    assert args.model == 'afrc'
    assert args.out == 'out'
    assert args.pdb_only is False
    assert args.seed is None


def test_long_and_short_options_agree():
    parser = build_ensemble_parser()
    short = parser.parse_args(['-s', 'ACDE', '-n', '7', '-m', 'afrc', '-o', 'x'])
    long = parser.parse_args(['--sequence', 'ACDE', '--number-of-conformers', '7', '--model', 'afrc', '--out', 'x'])
    assert vars(short) == vars(long)


def test_model_is_case_insensitive():
    assert build_ensemble_parser().parse_args(['-s', 'ACDE', '-m', 'AFRC']).model == 'afrc'
    assert build_ensemble_parser().parse_args(['-s', 'ACDE', '-m', 'WLC2']).model == 'wlc2'
    assert set(ENSEMBLE_MODELS) == {'afrc', 'fjc', 'frc', 'wlc', 'wlc2', 'saw', 'saw-nu'}


@pytest.mark.parametrize('argv', [
    [],                                          # no sequence
    ['-s', 'ACDE', '-n', '0'],                   # not positive
    ['-s', 'ACDE', '-n', '-5'],
    ['-s', 'ACDE', '-n', '2.5'],                 # not an integer
    ['-s', 'ACDE', '-n', 'many'],
    ['-s', 'ACDE', '-m', 'rouse'],               # unsupported model
    ['-s', 'ACDE', '--seed', 'abc'],
])
def test_bad_arguments_exit_with_status_2(argv, capsys):
    with pytest.raises(SystemExit) as error:
        ensemble_main(argv)
    assert error.value.code == 2
    assert 'afrc-ensemble' in capsys.readouterr().err


def test_help(capsys):
    with pytest.raises(SystemExit) as error:
        ensemble_main(['--help'])
    assert error.value.code == 0
    out = capsys.readouterr().out
    for option in ('--sequence', '--number-of-conformers', '--model', '--out', '--pdb-only', '--seed', 'examples:'):
        assert option in out


# ---------------------------------------------------------------------------
# a normal run
# ---------------------------------------------------------------------------
def test_writes_pdb_xtc_and_report(tmp_path, capsys):
    md = pytest.importorskip('mdtraj')
    assert _run(tmp_path, '-n', '300', '--seed', '5') == 0
    assert sorted(os.listdir(tmp_path)) == ['out.pdb', 'out.xtc', 'out_report.txt']

    traj = md.load(str(tmp_path / 'out.xtc'), top=str(tmp_path / 'out.pdb'))
    assert traj.n_frames == 300
    assert ''.join(residue.code for residue in traj.topology.residues) == P53

    # the ensemble is exactly the one the Python API gives for the same seed
    expected = afrc.AnalyticalFRC(P53).sample_conformations(300, seed=5)
    assert np.allclose(traj.xyz * 10, expected, atol=0.006)

    # the report is printed and saved, and says the ensemble is model-like
    printed = capsys.readouterr().out
    assert printed == _read(tmp_path / 'out_report.txt')
    assert 'Verdict: PASS' in printed
    assert 'Seed          : 5 ' in printed
    assert 'out.pdb, ' in printed and 'out.xtc' in printed
    assert P53 in printed


def test_pdb_only_writes_a_single_multimodel_file(tmp_path, capsys):
    md = pytest.importorskip('mdtraj')
    assert _run(tmp_path, '-n', '150', '--seed', '6', '--pdb-only') == 0
    assert sorted(os.listdir(tmp_path)) == ['out.pdb', 'out_report.txt']

    traj = md.load(str(tmp_path / 'out.pdb'))
    assert traj.n_frames == 150
    expected = afrc.AnalyticalFRC(P53).sample_conformations(150, seed=6)
    assert np.allclose(traj.xyz * 10, expected, atol=5e-4)
    assert 'Verdict: PASS' in capsys.readouterr().out


def test_pdb_only_does_not_need_mdtraj(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, 'mdtraj', None)
    monkeypatch.setitem(sys.modules, 'mdtraj.formats', None)
    assert _run(tmp_path, '-n', '120', '--pdb-only') == 0
    assert sorted(os.listdir(tmp_path)) == ['out.pdb', 'out_report.txt']


def test_seed_makes_runs_reproducible(tmp_path):
    first, second = tmp_path / 'a', tmp_path / 'b'
    first.mkdir()
    second.mkdir()
    for directory in (first, second):
        assert ensemble_main(['-s', P53, '-n', '120', '--seed', '42', '--pdb-only', '-o', str(directory / 'e')]) == 0
    assert _read(first / 'e.pdb') == _read(second / 'e.pdb')


def test_unseeded_runs_record_a_seed_that_reproduces_them(tmp_path):
    assert _run(tmp_path, '-n', '110', '--pdb-only') == 0
    seed = int(re.search(r'Seed\s+: (\d+)', _read(tmp_path / 'out_report.txt')).group(1))
    rerun = tmp_path / 'rerun'
    rerun.mkdir()
    assert ensemble_main(['-s', P53, '-n', '110', '--pdb-only', '--seed', str(seed), '-o', str(rerun / 'out')]) == 0
    assert _read(rerun / 'out.pdb') == _read(tmp_path / 'out.pdb')


def test_lowercase_sequence_and_extension_on_prefix(tmp_path):
    assert ensemble_main(['-s', P53.lower(), '-n', '110', '--pdb-only', '-o', str(tmp_path / 'named.pdb')]) == 0
    assert sorted(os.listdir(tmp_path)) == ['named.pdb', 'named_report.txt']


def test_short_run_is_reported_as_not_assessed(tmp_path, capsys):
    assert _run(tmp_path, '-n', '10', '--pdb-only') == 0
    assert 'Verdict: NOT ASSESSED' in capsys.readouterr().out


# ---------------------------------------------------------------------------
# errors: status 1, a clear message, and nothing written
# ---------------------------------------------------------------------------
def _assert_error(tmp_path, capsys, argv, message):
    assert ensemble_main(argv) == 1
    captured = capsys.readouterr()
    assert captured.err.startswith('afrc-ensemble: error: ')
    assert message in captured.err
    assert captured.out == ''
    assert os.listdir(tmp_path) == []


def test_invalid_sequence(tmp_path, capsys):
    _assert_error(tmp_path, capsys, ['-s', 'ACDXZ', '-o', str(tmp_path / 'out')], 'non-standard amino acids')


def test_missing_output_directory(tmp_path, capsys):
    _assert_error(tmp_path, capsys, ['-s', P53, '-o', str(tmp_path / 'nowhere' / 'out')], 'does not exist')


def test_missing_mdtraj_suggests_pdb_only(tmp_path, capsys, monkeypatch):
    monkeypatch.setitem(sys.modules, 'mdtraj', None)
    _assert_error(tmp_path, capsys, ['-s', P53, '-o', str(tmp_path / 'out')], '--pdb-only')


def test_too_many_conformations_for_a_multimodel_pdb(tmp_path, capsys):
    _assert_error(tmp_path, capsys, ['-s', P53, '-n', '10000', '--pdb-only', '-o', str(tmp_path / 'out')], 'at most 9999')


# ---------------------------------------------------------------------------
# installation
# ---------------------------------------------------------------------------
def test_console_script_is_declared():
    tomllib = pytest.importorskip('tomllib')
    pyproject = os.path.join(os.path.dirname(afrc.__file__), os.pardir, 'pyproject.toml')
    if not os.path.exists(pyproject):
        pytest.skip('not running from a source checkout')
    with open(pyproject, 'rb') as fh:
        scripts = tomllib.load(fh)['project']['scripts']
    assert scripts['afrc-ensemble'] == 'afrc.cli:ensemble_main'


def test_runs_as_a_module(tmp_path):
    result = subprocess.run([sys.executable, '-m', 'afrc.cli', '-s', 'ACDEFGHIKL', '-n', '120', '--pdb-only',
                             '--seed', '1', '-o', str(tmp_path / 'm')], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'Verdict: PASS' in result.stdout


# ---------------------------------------------------------------------------
# every polymer model
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('argv, model_line', [
    (['-m', 'afrc'], 'AFRC'),
    (['-m', 'fjc'], 'freely jointed chain (b = 3.8 A)'),
    (['-m', 'fjc', '--segment-length', '4.5'], 'freely jointed chain (b = 4.5 A)'),
    (['-m', 'frc', '--c-inf', '3'], 'freely rotating chain (b = 3.8 A, c_inf = 3)'),
    (['-m', 'wlc', '--lp', '4'], 'worm-like chain (lp = 4 A, aa_size = 3.8 A; Zhou)'),
    (['-m', 'wlc2', '--segment-length', '3.6'], "worm-like chain (lp = 3 A, aa_size = 3.6 A; O'Brien)"),
    (['-m', 'saw', '--prefactor', '6'], 'self-avoiding walk (prefactor = 6 A, nu = 0.598; Gaussian approximation)'),
    (['-m', 'saw-nu', '--nu', '0.45'], 'nu-dependent SAW (nu = 0.45, prefactor = 5.5 A; Gaussian approximation)'),
])
def test_every_model_runs_and_passes(tmp_path, capsys, argv, model_line):
    assert ensemble_main(['-s', P53, '-n', '200', '--seed', '3', '--pdb-only', '-o', str(tmp_path / 'out'), *argv]) == 0
    printed = capsys.readouterr().out
    assert f'Model         : {model_line}\n' in printed
    assert 'Verdict: PASS' in printed
    with open(tmp_path / 'out.pdb') as fh:
        assert model_line.split(' (')[0] in fh.readline()


def test_model_ensembles_match_the_python_api(tmp_path):
    from afrc.polymer_models.wlc import WormLikeChain

    assert ensemble_main(['-s', P53, '-n', '20', '--seed', '8', '--pdb-only', '-o', str(tmp_path / 'w'),
                          '-m', 'wlc', '--lp', '5']) == 0
    md = pytest.importorskip('mdtraj')
    expected = WormLikeChain(P53, lp=5.0).sample_conformations(20, seed=8)
    assert np.allclose(md.load(str(tmp_path / 'w.pdb')).xyz * 10, expected, atol=5e-4)


@pytest.mark.parametrize('argv, flag, users', [
    (['-m', 'afrc', '--lp', '4'], '--lp', 'wlc, wlc2'),
    (['-m', 'fjc', '--c-inf', '2'], '--c-inf', 'frc'),
    (['-m', 'saw', '--nu', '0.5'], '--nu', 'saw-nu'),
    (['-m', 'wlc', '--prefactor', '5'], '--prefactor', 'saw, saw-nu'),
    (['-m', 'afrc', '--segment-length', '3'], '--segment-length', 'fjc, frc, wlc, wlc2'),
])
def test_parameters_for_other_models_are_usage_errors(argv, flag, users, capsys):
    with pytest.raises(SystemExit) as error:
        ensemble_main(['-s', P53, *argv])
    assert error.value.code == 2
    assert f'{flag} does not apply' in capsys.readouterr().err
    # ...and the message says which models do use it
    with pytest.raises(SystemExit):
        ensemble_main(['-s', P53, *argv])
    assert f'(it is used by: {users})' in capsys.readouterr().err


@pytest.mark.parametrize('argv', [
    ['-m', 'wlc', '--lp', '0'],
    ['-m', 'wlc', '--lp', '-3'],
    ['-m', 'fjc', '--segment-length', 'long'],
    ['-m', 'saw-nu', '--nu', '1'],
    ['-m', 'saw-nu', '--nu', '0'],
    ['-m', 'saw', '--prefactor', '0'],
    ['-m', 'saw-nu', '--nu', 'half'],
])
def test_invalid_parameter_values_are_usage_errors(argv):
    with pytest.raises(SystemExit) as error:
        ensemble_main(['-s', P53, *argv])
    assert error.value.code == 2


def test_model_level_errors_are_reported_cleanly(tmp_path, capsys):
    """The O'Brien model needs the chain to be at least one persistence length long."""
    _assert_error(tmp_path, capsys, ['-s', 'AA', '-m', 'wlc2', '--lp', '20', '--pdb-only', '-o', str(tmp_path / 'out')],
                  'shorter than the persistence length')


def test_invalid_sequence_is_rejected_for_every_model(tmp_path, capsys):
    for model in ('fjc', 'saw-nu'):
        _assert_error(tmp_path, capsys, ['-s', 'ACDX', '-m', model, '--pdb-only', '-o', str(tmp_path / 'out')],
                      'non-standard amino acids')


def test_help_lists_every_model(capsys):
    with pytest.raises(SystemExit):
        ensemble_main(['--help'])
    out = capsys.readouterr().out
    for model in ('afrc', 'fjc', 'frc', 'wlc', 'wlc2', 'saw', 'saw-nu', '--segment-length', '--c-inf', '--lp', '--prefactor', '--nu'):
        assert model in out
