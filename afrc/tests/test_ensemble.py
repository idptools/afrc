"""
Tests for 3D ensemble generation (``AnalyticalFRC.sample_conformations()`` and
``save_ensemble()``) and the writers in ``afrc/ensemble.py``.

The statistical tests use fixed seeds, so they are deterministic; their
tolerances are several standard errors wide.
"""

import os
import sys

import numpy as np
import pytest
from scipy import stats

import afrc
from afrc.afrc import AFRCException
from afrc.config import AA_list
from afrc.ensemble import (THREE_LETTER_CODES, gaussian_chain_factor, sample_gaussian_chain,
                           write_ensemble, write_pdb, write_xtc)
from afrc.polymer import PolymerObject

P53 = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'


def _msd(p):
    return p.get_mean_squared_distance_map()


def _pairwise_distances(xyz):
    """[n_conformations x N x N] distance matrices."""
    return np.linalg.norm(xyz[:, :, None, :] - xyz[:, None, :, :], axis=-1)


# ---------------------------------------------------------------------------
# the mean-squared distance matrix
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('seed', range(3))
def test_mean_squared_distances_match_the_inter_residue_polymers(seed):
    seq = ''.join(np.random.default_rng(seed).choice(AA_list, 40))
    D = _msd(afrc.AnalyticalFRC(seq))
    assert np.array_equal(D, D.T)
    assert np.all(np.diag(D) == 0.0)
    for i in range(40):
        for j in range(i + 1, 40):
            assert D[i, j] == pytest.approx(PolymerObject(seq[i:j]).RMS_Re_scaling**2, rel=1e-12)


# ---------------------------------------------------------------------------
# sample_conformations: shape, seeding, validation
# ---------------------------------------------------------------------------
def test_sample_shape_and_type():
    xyz = afrc.AnalyticalFRC(P53).sample_conformations(n=25, seed=0)
    assert xyz.shape == (25, len(P53), 3)
    assert xyz.dtype == np.float64
    assert np.all(np.isfinite(xyz))


def test_conformations_are_centred():
    xyz = afrc.AnalyticalFRC(P53).sample_conformations(n=50, seed=0)
    assert np.allclose(xyz.mean(axis=1), 0.0, atol=1e-9)


def test_seeding_is_reproducible():
    p = afrc.AnalyticalFRC(P53)
    assert np.array_equal(p.sample_conformations(10, seed=3), p.sample_conformations(10, seed=3))
    assert np.array_equal(p.sample_conformations(10, seed=np.random.default_rng(3)), p.sample_conformations(10, seed=3))
    assert not np.array_equal(p.sample_conformations(10, seed=3), p.sample_conformations(10, seed=4))
    # no seed means a fresh draw each time
    assert not np.array_equal(p.sample_conformations(10), p.sample_conformations(10))


def test_seeding_does_not_touch_the_global_random_state():
    np.random.seed(12)
    expected = np.random.random()
    np.random.seed(12)
    afrc.AnalyticalFRC(P53).sample_conformations(10, seed=1)
    assert np.random.random() == expected


@pytest.mark.parametrize('bad', [0, -3, 2.5, '10', None, 'abc'])
def test_invalid_number_of_conformations_raises(bad):
    with pytest.raises(AFRCException):
        afrc.AnalyticalFRC(P53).sample_conformations(n=bad)


def test_integer_like_number_of_conformations_is_accepted():
    assert afrc.AnalyticalFRC(P53).sample_conformations(n=np.int64(4), seed=0).shape[0] == 4


def test_covariance_factor_is_built_once():
    p = afrc.AnalyticalFRC(P53)
    assert p._AnalyticalFRC__ensemble_factor is None
    p.sample_conformations(5, seed=0)
    factor = p._AnalyticalFRC__ensemble_factor
    p.sample_conformations(5, seed=1)
    assert p._AnalyticalFRC__ensemble_factor is factor


def test_ensemble_does_not_build_the_inter_residue_matrix():
    p = afrc.AnalyticalFRC(P53)
    p.sample_conformations(5, seed=0)
    assert p.matrix is False


def test_single_residue_is_a_bead_at_the_origin():
    xyz = afrc.AnalyticalFRC('W').sample_conformations(n=5, seed=0)
    assert xyz.shape == (5, 1, 3)
    assert np.all(xyz == 0.0)


@pytest.mark.parametrize('seq', ['AG' * 250, 'A' * 150 + 'G' + 'A' * 150,
                                 ''.join(('A' if b % 2 else 'G') * (b % 17 + 1) for b in range(40))])
def test_adversarial_compositions_are_still_realizable(seq):
    """Sequences that alternate or block the largest and smallest prefactors."""
    xyz = afrc.AnalyticalFRC(seq).sample_conformations(n=3, seed=0)
    assert np.all(np.isfinite(xyz))


# ---------------------------------------------------------------------------
# the ensemble reproduces the AFRC
# ---------------------------------------------------------------------------
@pytest.fixture(scope='module')
def p53_ensemble():
    p = afrc.AnalyticalFRC(P53)
    return p, p.sample_conformations(n=20000, seed=2024)


def test_mean_squared_distances_are_reproduced_for_every_pair(p53_ensemble):
    p, xyz = p53_ensemble
    msd = np.mean(_pairwise_distances(xyz[:5000])**2, axis=0)
    D = _msd(p)
    upper = np.triu_indices(len(P53), 1)
    # the relative standard error of <r^2> from 5000 samples is ~1.2%
    assert np.allclose(msd[upper], D[upper], rtol=0.06)
    assert np.mean(np.abs(msd[upper]/D[upper] - 1)) < 0.015


@pytest.mark.parametrize('i, j', [(0, 1), (10, 11), (0, 5), (7, 30), (0, 49)])
def test_distance_distributions_are_exactly_maxwell(p53_ensemble, i, j):
    """Each inter-residue distance follows the AFRC's Gaussian-chain (Maxwell) distribution."""
    p, xyz = p53_ensemble
    distances = np.linalg.norm(xyz[:, i] - xyz[:, j], axis=-1)
    scale = np.sqrt(_msd(p)[i, j] / 3.0)
    assert stats.kstest(distances, stats.maxwell(scale=scale).cdf).pvalue > 0.01
    # and the mean is the AFRC's mean inter-residue distance
    assert np.mean(distances) == pytest.approx(p.get_mean_interresidue_distance(i, j, 'distribution'), rel=0.02)


def test_end_to_end_convention(p53_ensemble):
    """The first-to-last distance is the (0, N-1) inter-residue distance, not the N-residue Re."""
    p, xyz = p53_ensemble
    n = len(P53)
    ree = np.mean(np.linalg.norm(xyz[:, -1] - xyz[:, 0], axis=-1))
    assert ree == pytest.approx(p.get_mean_interresidue_distance(0, n - 1, 'distribution'), rel=0.01)


def test_end_to_end_convention_is_visible_for_short_chains():
    p = afrc.AnalyticalFRC('ACDEFGHIKL')
    xyz = p.sample_conformations(n=40000, seed=5)
    ree = np.mean(np.linalg.norm(xyz[:, -1] - xyz[:, 0], axis=-1))
    assert ree == pytest.approx(p.get_mean_interresidue_distance(0, 9, 'distribution'), rel=0.01)
    assert ree < 0.97 * p.get_mean_end_to_end_distance('distribution')


def test_ensemble_is_isotropic(p53_ensemble):
    _, xyz = p53_ensemble
    end_to_end = xyz[:, -1] - xyz[:, 0]
    cov = np.cov(end_to_end.T)
    assert np.allclose(cov, np.eye(3) * np.trace(cov) / 3, atol=0.05 * np.trace(cov) / 3)


def test_adjacent_bead_spacing_is_the_afrc_neighbour_distribution(p53_ensemble):
    _, xyz = p53_ensemble
    spacing = np.linalg.norm(np.diff(xyz, axis=1), axis=-1)
    assert 5.5 < np.mean(spacing) < 6.2
    assert np.percentile(spacing, 5) < 3.0


@pytest.mark.parametrize('seq', [P53, 'G' * 60, 'A' * 20 + 'G' * 20 + 'P' * 20])
def test_radius_of_gyration_is_close_to_the_afrc(seq):
    p = afrc.AnalyticalFRC(seq)
    xyz = p.sample_conformations(n=5000, seed=9)
    rg = np.sqrt(np.mean(np.sum(xyz**2, axis=-1), axis=1))   # conformations are centred
    assert np.mean(rg) == pytest.approx(p.get_mean_radius_of_gyration(), rel=0.03)


def test_kirkwood_riseman_rh_is_reproduced(p53_ensemble):
    """1/<1/r_ij> over the ensemble is the AFRC's Kirkwood-Riseman Rh."""
    p, xyz = p53_ensemble
    d = _pairwise_distances(xyz[:4000])
    upper = np.triu_indices(len(P53), 1)
    rh = 1.0 / np.mean(1.0 / d[:, upper[0], upper[1]])
    assert rh == pytest.approx(p.get_mean_hydrodynamic_radius(), rel=0.01)


# ---------------------------------------------------------------------------
# gaussian_chain_factor / sample_gaussian_chain
# ---------------------------------------------------------------------------
def test_factor_reproduces_the_covariance_of_a_random_walk():
    n, b2 = 12, 14.44
    k = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    D = b2 * k
    A = gaussian_chain_factor(D)
    J = np.eye(n) - 1.0/n
    assert np.allclose(A @ A.T, -J @ D @ J / 6.0, atol=1e-10)


def test_factor_rejects_bad_matrices():
    with pytest.raises(AFRCException):
        gaussian_chain_factor(np.zeros((3, 4)))
    with pytest.raises(AFRCException):
        gaussian_chain_factor(np.array([[0.0, 1.0], [2.0, 0.0]]))
    # three points whose squared distances break the triangle inequality badly
    with pytest.raises(AFRCException):
        gaussian_chain_factor(np.array([[0.0, 1.0, 100.0], [1.0, 0.0, 1.0], [100.0, 1.0, 0.0]]))


def test_sample_gaussian_chain_shape():
    A = gaussian_chain_factor(14.44 * np.abs(np.subtract.outer(np.arange(6), np.arange(6))).astype(float))
    xyz = sample_gaussian_chain(A, 7, np.random.default_rng(0))
    assert xyz.shape == (7, 6, 3)


# ---------------------------------------------------------------------------
# PDB writer
# ---------------------------------------------------------------------------
def _read(path):
    with open(path) as fh:
        return fh.read().splitlines()


def test_pdb_format(tmp_path):
    xyz = afrc.AnalyticalFRC(P53).sample_conformations(1, seed=0)[0]
    path = str(tmp_path / 'one.pdb')
    write_pdb(xyz, P53, path, remark='a remark that is long enough that it has to be wrapped onto a second REMARK line to fit')
    lines = _read(path)

    assert all(len(line) <= 80 for line in lines)
    remarks = [line for line in lines if line.startswith('REMARK')]
    assert len(remarks) == 2

    atoms = [line for line in lines if line.startswith('ATOM')]
    assert len(atoms) == len(P53)
    for index, (line, residue) in enumerate(zip(atoms, P53), start=1):
        assert int(line[6:11]) == index
        assert line[12:16] == ' CA '
        assert line[17:20] == THREE_LETTER_CODES[residue]
        assert line[21] == 'A'
        assert int(line[22:26]) == index
        assert line[76:78] == ' C'
        coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        assert np.allclose(coords, xyz[index - 1], atol=5e-4)

    conect = [line for line in lines if line.startswith('CONECT')]
    assert len(conect) == len(P53) - 1
    assert conect[0] == 'CONECT    1    2'
    assert any(line.startswith('TER') for line in lines)
    assert lines[-1] == 'END'


def test_pdb_without_remark_or_bonds(tmp_path):
    path = str(tmp_path / 'bead.pdb')
    write_pdb(np.zeros((1, 3)), 'G', path)
    lines = _read(path)
    assert not any(line.startswith(('REMARK', 'CONECT')) for line in lines)
    assert sum(line.startswith('ATOM') for line in lines) == 1


@pytest.mark.parametrize('coordinates, sequence', [
    (np.zeros((2, 5, 3)), 'AAAAA'),       # more than one conformation
    (np.zeros((4, 3)), 'AAAAA'),          # wrong number of residues
    (np.zeros((5, 2)), 'AAAAA'),          # not 3D
    (np.zeros((5, 3)), 'AAXAA'),          # non-standard residue
    (np.zeros((0, 3)), ''),               # empty
    (np.full((2, 3), 1000.0), 'AA'),      # too large for the PDB columns
])
def test_pdb_writer_rejects_bad_input(tmp_path, coordinates, sequence):
    with pytest.raises(AFRCException):
        write_pdb(coordinates, sequence, str(tmp_path / 'bad.pdb'))


def test_pdb_writer_rejects_too_many_residues(tmp_path):
    with pytest.raises(AFRCException):
        write_pdb(np.zeros((10000, 3)), 'A' * 10000, str(tmp_path / 'big.pdb'))


# ---------------------------------------------------------------------------
# XTC writer and the PDB/XTC pair
# ---------------------------------------------------------------------------
def test_xtc_round_trip(tmp_path):
    md = pytest.importorskip('mdtraj')
    xyz = afrc.AnalyticalFRC(P53).sample_conformations(30, seed=1)
    path = str(tmp_path / 'traj.xtc')
    write_xtc(xyz, path)
    with md.formats.XTCTrajectoryFile(path) as fh:
        read, _, _, _ = fh.read()
    assert read.shape == xyz.shape
    # XTC stores 0.001 nm = 0.01 A
    assert np.allclose(read * 10, xyz, atol=0.006)


def test_xtc_writer_rejects_bad_shapes(tmp_path):
    pytest.importorskip('mdtraj')
    with pytest.raises(AFRCException):
        write_xtc(np.zeros((5, 3)), str(tmp_path / 'bad.xtc'))


def test_save_ensemble_writes_a_loadable_pair(tmp_path):
    md = pytest.importorskip('mdtraj')
    p = afrc.AnalyticalFRC(P53)
    prefix = str(tmp_path / 'p53')
    xyz = p.save_ensemble(prefix, n=40, seed=11)

    assert np.array_equal(xyz, p.sample_conformations(40, seed=11))
    traj = md.load(prefix + '.xtc', top=prefix + '.pdb')
    assert traj.n_frames == 40
    assert ''.join(residue.code for residue in traj.topology.residues) == P53
    assert traj.topology.n_bonds == len(P53) - 1
    assert np.allclose(traj.xyz * 10, xyz, atol=0.006)

    pdb = md.load(prefix + '.pdb')
    assert np.allclose(pdb.xyz[0] * 10, xyz[0], atol=5e-4)
    assert afrc.__version__ in ' '.join(_read(prefix + '.pdb')[:2])


@pytest.mark.parametrize('name', ['ens.pdb', 'ens.xtc', 'ens.PDB', 'ens'])
def test_save_ensemble_file_names(tmp_path, name):
    pytest.importorskip('mdtraj')
    afrc.AnalyticalFRC('ACDEFG').save_ensemble(str(tmp_path / name), n=2, seed=0)
    assert sorted(os.listdir(tmp_path)) == ['ens.pdb', 'ens.xtc']


def test_missing_mdtraj_fails_cleanly(tmp_path, monkeypatch):
    """Without mdtraj nothing is written and the error says how to fix it."""
    monkeypatch.setitem(sys.modules, 'mdtraj', None)
    monkeypatch.setitem(sys.modules, 'mdtraj.formats', None)
    with pytest.raises(AFRCException, match='mdtraj'):
        afrc.AnalyticalFRC('ACDEFG').save_ensemble(str(tmp_path / 'ens'), n=2, seed=0)
    assert os.listdir(tmp_path) == []


def test_write_ensemble_checks_the_pdb_before_writing_anything(tmp_path):
    pytest.importorskip('mdtraj')
    xyz = np.zeros((3, 4, 3))
    xyz[0, 0, 0] = 5000.0
    with pytest.raises(AFRCException):
        write_ensemble(xyz, 'ACDE', str(tmp_path / 'e.pdb'), str(tmp_path / 'e.xtc'))
    assert os.listdir(tmp_path) == []


# ---------------------------------------------------------------------------
# multi-model PDB writer
# ---------------------------------------------------------------------------
def test_multimodel_pdb_format(tmp_path):
    from afrc.ensemble import write_multimodel_pdb

    xyz = afrc.AnalyticalFRC(P53).sample_conformations(7, seed=2)
    path = str(tmp_path / 'multi.pdb')
    write_multimodel_pdb(xyz, P53, path, remark='seven conformations')
    lines = _read(path)

    assert all(len(line) <= 80 for line in lines)
    assert lines[0] == 'REMARK   1 seven conformations'
    models = [line for line in lines if line.startswith('MODEL')]
    assert models == [f'MODEL     {i:4d}' for i in range(1, 8)]
    assert sum(line == 'ENDMDL' for line in lines) == 7
    assert sum(line.startswith('ATOM') for line in lines) == 7 * len(P53)
    # bonds are written once, after the last model
    conect = [i for i, line in enumerate(lines) if line.startswith('CONECT')]
    assert len(conect) == len(P53) - 1
    assert min(conect) > max(i for i, line in enumerate(lines) if line == 'ENDMDL')
    assert lines[-1] == 'END'


def test_multimodel_pdb_round_trip(tmp_path):
    md = pytest.importorskip('mdtraj')
    from afrc.ensemble import write_multimodel_pdb

    xyz = afrc.AnalyticalFRC(P53).sample_conformations(12, seed=3)
    path = str(tmp_path / 'multi.pdb')
    write_multimodel_pdb(xyz, P53, path)
    traj = md.load(path)
    assert traj.n_frames == 12
    assert traj.topology.n_bonds == len(P53) - 1
    assert ''.join(residue.code for residue in traj.topology.residues) == P53
    assert np.allclose(traj.xyz * 10, xyz, atol=5e-4)


def test_multimodel_pdb_limits(tmp_path):
    from afrc.ensemble import write_multimodel_pdb

    with pytest.raises(AFRCException, match='9999 models'):
        write_multimodel_pdb(np.zeros((10000, 2, 3)), 'AG', str(tmp_path / 'many.pdb'))
    xyz = np.zeros((3, 2, 3))
    xyz[2, 1, 0] = 2000.0   # only the last model is too large
    with pytest.raises(AFRCException, match='XTC'):
        write_multimodel_pdb(xyz, 'AG', str(tmp_path / 'far.pdb'))
    assert os.listdir(tmp_path) == []


def test_save_ensemble_pdb_only(tmp_path, monkeypatch):
    """pdb_only writes one multi-model file and needs no mdtraj."""
    monkeypatch.setitem(sys.modules, 'mdtraj', None)
    monkeypatch.setitem(sys.modules, 'mdtraj.formats', None)
    p = afrc.AnalyticalFRC('ACDEFGHIKL')
    xyz = p.save_ensemble(str(tmp_path / 'ens.pdb'), n=5, seed=4, pdb_only=True)
    assert os.listdir(tmp_path) == ['ens.pdb']
    assert np.array_equal(xyz, p.sample_conformations(5, seed=4))
    assert sum(line.startswith('MODEL') for line in _read(str(tmp_path / 'ens.pdb'))) == 5


@pytest.mark.parametrize('sequence, n_models, remark', [
    ('W', 3, ''),
    ('GS', 4, 'short'),
    ('ACDEFGHIKLMNPQRSTVWY', 25, 'a remark long enough that it has to be wrapped onto a second REMARK line to fit'),
    (P53, 60, 'x'),
])
def test_multimodel_pdb_text_matches_line_by_line_formatting(sequence, n_models, remark):
    """
    The multi-model writer fills in one pre-built block per model rather than
    formatting every line; the text must be exactly what line-by-line formatting
    (write_pdb's layout) gives, including rounding and the column limits.
    """
    from afrc.ensemble import _atom_lines, _conect_lines, _multimodel_pdb_chunks, _remark_lines

    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(n_models, len(sequence), 3)) * rng.choice([0.001, 1.0, 50.0, 400.0], size=(n_models, len(sequence), 3))
    xyz[0, 0] = [-0.0004, 999.999, -999.999]

    lines = _remark_lines(remark)
    for model, frame in enumerate(xyz, start=1):
        lines += [f'MODEL     {model:4d}'] + _atom_lines(frame, sequence) + ['ENDMDL']
    expected = '\n'.join(lines + _conect_lines(sequence) + ['END']) + '\n'

    assert ''.join(_multimodel_pdb_chunks(xyz, sequence, remark)) == expected
