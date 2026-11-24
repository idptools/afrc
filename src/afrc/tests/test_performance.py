"""
Tests for the performance-related behaviour of ``AnalyticalFRC`` and ``PolymerObject``.

These pin three things: (1) single-pair queries do not build the [n x n]
inter-residue matrix, (2) matrix entries do not cache their distributions (which
used to cost gigabytes of memory), and (3) none of this changes any result.
"""

import tracemalloc

import numpy as np
import pytest

import afrc
from afrc.config import AA_list, RIJ_R0, RIJ_RMS_R0, RG_R0, RG_X0
from afrc.polymer import PolymerObject


# ---------------------------------------------------------------------------
# PolymerObject caching
# ---------------------------------------------------------------------------
def test_caching_is_on_by_default(all_aa):
    po = PolymerObject(all_aa)
    assert po.cache_distributions is True
    assert po.get_end_to_end_distribution()[1] is po.get_end_to_end_distribution()[1]
    assert po.get_radius_of_gyration_distribution()[1] is po.get_radius_of_gyration_distribution()[1]


def test_uncached_polymer_object_stores_nothing(all_aa):
    po = PolymerObject(all_aa, cache_distributions=False)
    first_re = po.get_end_to_end_distribution()
    first_rg = po.get_radius_of_gyration_distribution()
    po.get_mean_end_to_end_distance('distribution')
    po.get_mean_radius_of_gyration('distribution')
    po.sample_end_to_end_distribution(10)
    po.sample_radius_of_gyration_distribution(10)

    # nothing was kept...
    assert po._PolymerObject__p_of_Re_R is False
    assert po._PolymerObject__p_of_Rg_R is False

    # ...so each call builds fresh (but identical) arrays
    second_re = po.get_end_to_end_distribution()
    assert second_re[1] is not first_re[1]
    assert np.array_equal(second_re[0], first_re[0]) and np.array_equal(second_re[1], first_re[1])
    second_rg = po.get_radius_of_gyration_distribution()
    assert np.array_equal(second_rg[1], first_rg[1])


@pytest.mark.parametrize('seq', ['', 'A', 'GS', 'ACDEFGHIKLMNPQRSTVWY', 'MEEPQSDPSVEPPLSQETFSDLWKLLPENN'])
def test_caching_does_not_change_results(seq):
    cached = PolymerObject(seq)
    uncached = PolymerObject(seq, cache_distributions=False)

    for getter in ('get_end_to_end_distribution', 'get_radius_of_gyration_distribution'):
        for a, b in zip(getattr(cached, getter)(), getattr(uncached, getter)()):
            assert np.array_equal(a, b)

    for mode in ('scaling law', 'distribution'):
        assert cached.get_mean_end_to_end_distance(mode) == uncached.get_mean_end_to_end_distance(mode)
        assert cached.get_mean_radius_of_gyration(mode) == uncached.get_mean_radius_of_gyration(mode)

    np.random.seed(4)
    a = (cached.sample_end_to_end_distribution(50), cached.sample_radius_of_gyration_distribution(50))
    np.random.seed(4)
    b = (uncached.sample_end_to_end_distribution(50), uncached.sample_radius_of_gyration_distribution(50))
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


# ---------------------------------------------------------------------------
# composition-weighted prefactors (counted once per residue type)
# ---------------------------------------------------------------------------
def _reference_prefactors(seq):
    """The original implementation: every residue type counted, for every prefactor."""
    n = len(seq)
    r0_rms = r0 = rg_r0 = x0 = 0
    for aa in AA_list:
        r0_rms = r0_rms + (seq.count(aa)/float(n))*RIJ_RMS_R0[aa]
        r0 = r0 + (seq.count(aa)/float(n))*RIJ_R0[aa]
        rg_r0 = rg_r0 + (seq.count(aa)/float(n))*RG_R0[aa]
        x0 = x0 + (seq.count(aa)/float(n))*(RG_X0[aa] + 0.005)
    return r0_rms, r0, rg_r0, x0


@pytest.mark.parametrize('seed', range(6))
def test_prefactors_are_bitwise_identical_to_reference(seed):
    """Skipping absent residues only drops exact +0.0 terms, so nothing may change."""
    rng = np.random.default_rng(seed)
    # vary both length and how many residue types are present
    alphabet = rng.choice(AA_list, size=rng.integers(1, 21), replace=False)
    seq = ''.join(rng.choice(alphabet, size=rng.integers(1, 300)))

    po = PolymerObject(seq)
    r0_rms, r0, rg_r0, x0 = _reference_prefactors(seq)
    assert po._PolymerObject__R0_RMS == r0_rms
    assert po._PolymerObject__R0 == r0
    assert po._PolymerObject__RG_R0 == rg_r0
    assert po._PolymerObject__X0 == x0


# ---------------------------------------------------------------------------
# single-pair queries never build the whole matrix
# ---------------------------------------------------------------------------
SINGLE_PAIR_CALLS = {
    'distance distribution': lambda p: p.get_interresidue_distance_distribution(3, 17),
    'mean distance (scaling law)': lambda p: p.get_mean_interresidue_distance(3, 17),
    'mean distance (distribution)': lambda p: p.get_mean_interresidue_distance(3, 17, 'distribution'),
    'segment Rg (scaling law)': lambda p: p.get_mean_interresidue_radius_of_gyration(3, 17),
    'segment Rg (distribution)': lambda p: p.get_mean_interresidue_radius_of_gyration(3, 17, 'distribution'),
    'sampling': lambda p: p.sample_inter_residue_distance_distribution(3, 17, n=20),
    'contact fraction': lambda p: p.get_contact_fraction(3, 17, 12.0),
    'PRE profile': lambda p: p.get_pre_profile(3, sample_size=50),
}


@pytest.mark.parametrize('name', SINGLE_PAIR_CALLS)
def test_single_pair_queries_do_not_build_the_matrix(name, test_seq):
    p = afrc.AnalyticalFRC(test_seq)
    SINGLE_PAIR_CALLS[name](p)
    assert p.matrix is False


@pytest.mark.parametrize('name', [n for n in SINGLE_PAIR_CALLS if n != 'sampling'])
def test_single_pair_results_do_not_depend_on_the_matrix(name, test_seq):
    """Before and after the matrix exists, every single-pair answer is bitwise identical."""
    p = afrc.AnalyticalFRC(test_seq)
    np.random.seed(9)
    before = SINGLE_PAIR_CALLS[name](p)
    p.get_distance_map()
    assert p.matrix is not False
    np.random.seed(9)
    after = SINGLE_PAIR_CALLS[name](p)

    flat_before, flat_after = _flatten(before), _flatten(after)
    assert len(flat_before) == len(flat_after)
    for a, b in zip(flat_before, flat_after):
        assert np.array_equal(a, b)


def _flatten(x):
    """Flatten nested lists/tuples of arrays and scalars into a list of float arrays."""
    if isinstance(x, (list, tuple)):
        out = []
        for item in x:
            out.extend(_flatten(item))
        return out
    return [np.asarray(x, dtype=float)]


def test_single_pair_segment_matches_the_matrix_entry(test_seq):
    p = afrc.AnalyticalFRC(test_seq)
    fresh = p._AnalyticalFRC__get_pair(40, 12)
    assert fresh.nres == 28
    assert fresh.cache_distributions is False
    p.get_distance_map()
    from_matrix = p._AnalyticalFRC__get_pair(40, 12)
    assert from_matrix is p.matrix[12][40]
    assert fresh.get_mean_end_to_end_distance() == from_matrix.get_mean_end_to_end_distance()
    assert np.array_equal(fresh.get_end_to_end_distribution()[1], from_matrix.get_end_to_end_distribution()[1])


# ---------------------------------------------------------------------------
# the matrix does not cache distributions
# ---------------------------------------------------------------------------
def test_matrix_entries_do_not_cache_distributions(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    p.get_contact_map(10.0)
    p.get_distance_map('distribution')
    p.get_internal_scaling('distribution')
    n = len(all_aa)
    for i in range(n):
        for j in range(n):
            entry = p.matrix[i][j]
            assert entry.cache_distributions is False
            assert entry._PolymerObject__p_of_Re_R is False


def test_full_chain_polymer_still_caches(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    assert p.full_seq_PO.cache_distributions is True
    assert p.get_end_to_end_distribution()[1] is p.get_end_to_end_distribution()[1]
    assert p.get_radius_of_gyration_distribution()[1] is p.get_radius_of_gyration_distribution()[1]


def test_bulk_maps_do_not_retain_memory(test_seq):
    """
    Regression: every pair's P(r) used to be cached, so a contact map or a
    distribution-mode distance map left ~100 MB allocated at 80 residues (and
    6.8 GB at 400). Nothing should now be retained beyond the matrix itself.
    """
    p = afrc.AnalyticalFRC(test_seq[:80])
    p.get_distance_map()   # build the matrix first; we only measure the maps

    tracemalloc.start()
    try:
        p.get_contact_map(10.0)
        p.get_distance_map('distribution')
        retained, _ = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert retained < 10e6


def test_repeated_bulk_calls_give_identical_results(all_aa):
    p = afrc.AnalyticalFRC(all_aa)
    assert np.array_equal(p.get_contact_map(8.0), p.get_contact_map(8.0))
    assert np.array_equal(p.get_distance_map('distribution'), p.get_distance_map('distribution'))
