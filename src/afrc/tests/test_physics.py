"""
Physical-consistency tests: invariances of the AFRC, its limiting behaviour, and
agreement between the models where theory says they should agree. Several of
these pin statements made in the documentation.
"""

import numpy as np
import pytest

import afrc
from afrc.config import AA_list, RIJ_R0, RIJ_RMS_R0
from afrc.polymer import PolymerObject
from afrc.polymer_models.fjc import FreelyJointedChain
from afrc.polymer_models.frc import FreelyRotatingChain
from afrc.polymer_models.nudep_saw import NuDepSAW
from afrc.polymer_models.saw import SAW
from afrc.polymer_models.wlc import WormLikeChain
from afrc.polymer_models.wlc2 import WormLikeChain2


def _random_sequence(seed, n):
    rng = np.random.default_rng(seed)
    return ''.join(rng.choice(AA_list, size=n))


# ---------------------------------------------------------------------------
# invariances of the AFRC
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('seed', range(3))
def test_whole_chain_properties_depend_only_on_composition(seed):
    """Shuffling a sequence leaves every whole-chain quantity bitwise unchanged."""
    seq = _random_sequence(seed, 90)
    shuffled = ''.join(np.random.default_rng(seed + 100).permutation(list(seq)))
    assert shuffled != seq

    a, b = afrc.AnalyticalFRC(seq), afrc.AnalyticalFRC(shuffled)
    for mode in ('scaling law', 'distribution'):
        assert a.get_mean_end_to_end_distance(mode) == b.get_mean_end_to_end_distance(mode)
        assert a.get_mean_radius_of_gyration(mode) == b.get_mean_radius_of_gyration(mode)
    assert a.get_mean_hydrodynamic_radius('nygaard') == b.get_mean_hydrodynamic_radius('nygaard')
    for x, y in zip(a.get_end_to_end_distribution(), b.get_end_to_end_distribution()):
        assert np.array_equal(x, y)


def test_sequence_order_matters_for_inter_residue_quantities():
    """...but not for inter-residue ones, which see the local composition."""
    blocky = afrc.AnalyticalFRC('A' * 20 + 'G' * 20)
    mixed = afrc.AnalyticalFRC('AG' * 20)
    assert blocky.get_mean_end_to_end_distance() == mixed.get_mean_end_to_end_distance()
    assert blocky.get_mean_interresidue_distance(0, 10) != mixed.get_mean_interresidue_distance(0, 10)
    assert blocky.get_mean_hydrodynamic_radius() != mixed.get_mean_hydrodynamic_radius()


def test_homopolymer_maps_depend_only_on_separation():
    """For a homopolymer every pair at the same |i - j| sees the same segment."""
    n = 25
    p = afrc.AnalyticalFRC('S' * n)
    dm = p.get_distance_map()
    cm = p.get_contact_map(9.0)
    for k in range(1, n):
        assert np.all(np.diagonal(dm, offset=k) == dm[0, k])
        assert np.all(np.diagonal(cm, offset=k) == cm[0, k])


def test_homopolymer_internal_scaling_is_the_scaling_law():
    n = 40
    scaling = afrc.AnalyticalFRC('T' * n).get_internal_scaling()
    k = np.arange(1, n)
    assert np.array_equal(scaling[:, 0], k)
    assert np.allclose(scaling[:, 1], RIJ_R0['T'] * np.sqrt(k), rtol=1e-12)


@pytest.mark.parametrize('seed', range(3))
def test_heteropolymer_internal_scaling_exponent_is_one_half(seed):
    scaling = afrc.AnalyticalFRC(_random_sequence(seed, 120)).get_internal_scaling()
    slope = np.polyfit(np.log(scaling[:, 0]), np.log(scaling[:, 1]), 1)[0]
    assert slope == pytest.approx(0.5, abs=0.005)


@pytest.mark.parametrize('seed', range(3))
def test_pair_distribution_is_the_segment_distribution(seed):
    seq = _random_sequence(seed, 60)
    p = afrc.AnalyticalFRC(seq)
    i, j = 7, 41
    expected = PolymerObject(seq[i:j]).get_end_to_end_distribution()
    for x, y in zip(p.get_interresidue_distance_distribution(j, i), expected):
        assert np.array_equal(x, y)


# ---------------------------------------------------------------------------
# size ordering and scaling of the AFRC
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('seed', range(3))
@pytest.mark.parametrize('n', [10, 40, 150])
def test_rh_rg_re_ordering(seed, n):
    """For all but the shortest chains Rh < Rg < Re."""
    p = afrc.AnalyticalFRC(_random_sequence(seed, n))
    assert p.get_mean_hydrodynamic_radius() < p.get_mean_radius_of_gyration() < p.get_mean_end_to_end_distance()


def test_rg_over_rh_approaches_the_kirkwood_limit():
    """
    For an ideal chain the Kirkwood-Riseman ratio sqrt(<Rg^2>)/Rh tends to
    8/(3 sqrt(pi)) ~ 1.505 from below as the chain grows.
    """
    limit = 8 / (3 * np.sqrt(np.pi))
    ratios = []
    for n in (25, 50, 100, 200):
        p = afrc.AnalyticalFRC('G' * n)
        rms_rg = p.full_seq_PO.RMS_Re_scaling / np.sqrt(6)
        ratios.append(rms_rg / p.get_mean_hydrodynamic_radius())
    assert np.all(np.diff(ratios) > 0)
    assert np.all(np.asarray(ratios) < limit)
    assert ratios[-1] > 0.9 * limit


def test_rh_grows_just_under_n_to_the_half():
    """Rh approaches N^0.5 from below, so its local exponent sits just under 0.5."""
    rh = [afrc.AnalyticalFRC('G' * n).get_mean_hydrodynamic_radius() for n in (100, 200)]
    exponent = np.log(rh[1] / rh[0]) / np.log(2)
    assert 0.4 < exponent < 0.5


def test_docs_rg_over_rh_for_typical_idrs(test_seq):
    """The docs quote Rg/Rh ~ 1.3-1.4 for typical IDR lengths."""
    p = afrc.AnalyticalFRC(test_seq)
    assert 1.3 < p.get_mean_radius_of_gyration() / p.get_mean_hydrodynamic_radius() < 1.4


# ---------------------------------------------------------------------------
# agreement between models
# ---------------------------------------------------------------------------
@pytest.mark.parametrize('seq', ['MASNDYTQQATQSYGAYPTQPGQGYSQQSS', 'G' * 80, 'KE' * 77])
def test_nu_dependent_saw_at_one_half_reproduces_the_afrc(seq):
    """With prefactor = R0_rms, the SAW-nu at nu = 0.5 has exactly the AFRC's RMS size."""
    p = afrc.AnalyticalFRC(seq)
    prefactor = p.full_seq_PO.RMS_Re_scaling / np.sqrt(len(seq))
    model = NuDepSAW(seq)
    assert model.get_root_mean_squared_end_to_end_distance(0.5, prefactor) == pytest.approx(p.full_seq_PO.RMS_Re_scaling, rel=1e-8)
    assert model.get_mean_end_to_end_distance(0.5, prefactor) == pytest.approx(p.get_mean_end_to_end_distance('distribution'), rel=0.01)


@pytest.mark.parametrize('n', [50, 150, 300])
def test_afrc_matches_a_worm_like_chain_with_lp_about_5(n):
    """The docs say matching the AFRC needs lp ~ 5 A; the ideal-chain value is R0_rms^2/(2b)."""
    lp = RIJ_RMS_R0['G']**2 / (2 * 3.8)
    assert 4.5 < lp < 5.5
    wlc = WormLikeChain2('G' * n, lp=lp)
    assert PolymerObject('G' * n).RMS_Re_scaling == pytest.approx(wlc.get_root_mean_squared_end_to_end_distance(), rel=0.01)


@pytest.mark.parametrize('n', [50, 154])
def test_docs_comparison_between_afrc_and_worm_like_chain(n, test_seq):
    """At lp = 3 A the AFRC is ~30% larger than the WLC; at lp = 3-4 A the WLC is 10-25% smaller."""
    seq = (test_seq * 2)[:n]
    afrc_mean = afrc.AnalyticalFRC(seq).get_mean_end_to_end_distance('distribution')
    wlc3 = WormLikeChain(seq, lp=3.0).get_mean_end_to_end_distance()
    wlc4 = WormLikeChain(seq, lp=4.0).get_mean_end_to_end_distance()
    assert 1.25 < afrc_mean / wlc3 < 1.40
    assert 0.75 < wlc3 / afrc_mean < 0.90
    assert 0.75 < wlc4 / afrc_mean < 0.90


@pytest.mark.parametrize('n', [30, 60, 120, 300])
def test_zhou_and_obrien_worm_like_chains_agree(n):
    zhou, obrien = WormLikeChain('A' * n), WormLikeChain2('A' * n)
    assert zhou.get_mean_end_to_end_distance() == pytest.approx(obrien.get_mean_end_to_end_distance(), rel=0.01)
    assert zhou.get_root_mean_squared_end_to_end_distance() == pytest.approx(obrien.get_root_mean_squared_end_to_end_distance(), rel=0.005)


def test_freely_rotating_chain_with_c_inf_one_converges_to_the_freely_jointed_chain():
    """Both are ideal chains of N bonds of length b; they differ only by finite extensibility."""
    gaps = []
    for n in (20, 100, 400):
        frc = FreelyRotatingChain('A' * n, c_inf=1.0).get_root_mean_squared_end_to_end_distance()
        fjc = FreelyJointedChain('A' * n).get_root_mean_squared_end_to_end_distance()
        assert frc > fjc
        gaps.append(frc / fjc - 1)
    assert np.all(np.diff(gaps) < 0)
    assert gaps[-1] < 0.005


def test_excluded_volume_models_order_by_nu(test_seq):
    """Collapsed < ideal < good solvent, and the fixed SAW sits at the top."""
    model = NuDepSAW(test_seq)
    sizes = [model.get_mean_end_to_end_distance(nu=nu) for nu in (0.35, 0.45, 0.5, 0.55, 0.598)]
    assert np.all(np.diff(sizes) > 0)
    assert SAW(test_seq).get_mean_end_to_end_distance() == pytest.approx(sizes[-1], rel=0.01)
    assert SAW(test_seq).get_mean_end_to_end_distance() > afrc.AnalyticalFRC(test_seq).get_mean_end_to_end_distance()


@pytest.mark.parametrize('n', [50, 200])
def test_saw_rg_to_re_ratio_is_universal(n):
    """The SAW's Rg/Re ratio does not depend on chain length or prefactor."""
    model = SAW('A' * n)
    ratio = model.get_mean_radius_of_gyration() / model.get_root_mean_squared_end_to_end_distance()
    ratio_other = model.get_mean_radius_of_gyration(prefactor=4.0) / model.get_root_mean_squared_end_to_end_distance(prefactor=4.0)
    assert ratio == pytest.approx(ratio_other, rel=1e-12)
    assert ratio == pytest.approx(SAW('A' * 10).get_mean_radius_of_gyration() / SAW('A' * 10).get_root_mean_squared_end_to_end_distance(), rel=1e-12)


def test_stiffer_chains_are_larger():
    lengths = [WormLikeChain2('A' * 80, lp=lp).get_mean_end_to_end_distance() for lp in (2.0, 3.0, 5.0, 10.0)]
    assert np.all(np.diff(lengths) > 0)
    lengths = [FreelyRotatingChain('A' * 80, c_inf=c).get_mean_end_to_end_distance() for c in (1.0, 2.0, 4.0, 9.0)]
    assert np.all(np.diff(lengths) > 0)
