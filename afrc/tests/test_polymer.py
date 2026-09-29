"""
Tests for ``afrc/polymer.py`` - the internal ``PolymerObject`` class.
"""

import numpy as np
import pytest

from afrc.exceptions import AFRCException
from afrc.polymer import PolymerObject


def test_zero_length_polymer_object():
    po = PolymerObject('')
    assert po.zero_length is True
    assert np.all(po.sample_end_to_end_distribution(10) == 0.0)
    assert np.all(po.sample_radius_of_gyration_distribution(10) == 0.0)


def test_zero_length_distributions_are_degenerate():
    """Regression: these previously crashed instead of returning a degenerate PMF."""
    po = PolymerObject('')
    for dist, prob in (po.get_end_to_end_distribution(),
                       po.get_radius_of_gyration_distribution()):
        assert dist == pytest.approx([0.0])
        assert prob == pytest.approx([1.0])


def test_apparent_rms_bond_length_positive(all_aa):
    """Regression test: this method previously crashed (wrong attribute name)."""
    po = PolymerObject(all_aa)
    b = po.compute_apparent_rms_bond_length()
    assert b > 0
    assert np.isfinite(b)


def test_apparent_rms_bond_length_needs_two_residues():
    """A one-residue chain has no bonds, so this is undefined rather than inf."""
    with pytest.raises(AFRCException):
        PolymerObject('A').compute_apparent_rms_bond_length()


def test_end_to_end_distribution_is_valid_pmf(all_aa):
    po = PolymerObject(all_aa)
    dist, prob = po.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_rg_distribution_is_valid_pmf(all_aa):
    po = PolymerObject(all_aa)
    dist, prob = po.get_radius_of_gyration_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_end_to_end_scaling_vs_distribution(all_aa):
    po = PolymerObject(all_aa)
    re_law = po.get_mean_end_to_end_distance('scaling law')
    re_dist = po.get_mean_end_to_end_distance('distribution')
    assert re_law > 0
    assert re_dist == pytest.approx(re_law, rel=0.05)


def test_mean_radius_of_gyration_scaling_law_matches_distribution(all_aa):
    """The calibrated Rg prefactor should reproduce the mean of the Rg distribution."""
    po = PolymerObject(all_aa)
    rg_law = po.get_mean_radius_of_gyration('scaling law')
    rg_dist = po.get_mean_radius_of_gyration('distribution')
    assert rg_law == pytest.approx(rg_dist, rel=0.01)


def test_invalid_calculation_mode_raises(all_aa):
    """PolymerObject should raise rather than silently returning None."""
    po = PolymerObject(all_aa)
    with pytest.raises(AFRCException):
        po.get_mean_end_to_end_distance('not-a-mode')
    with pytest.raises(AFRCException):
        po.get_mean_radius_of_gyration('not-a-mode')


def test_mean_radius_of_gyration_distribution_positive(all_aa):
    po = PolymerObject(all_aa)
    assert po.get_mean_radius_of_gyration('distribution') > 0


def test_sampling_sizes(all_aa):
    po = PolymerObject(all_aa)
    assert len(po.sample_end_to_end_distribution(64)) == 64
    assert len(po.sample_radius_of_gyration_distribution(64)) == 64


def test_mean_inverse_distance_matches_numerical_integration(all_aa):
    """<1/r> = sqrt(6 / (pi <r^2>)) should agree with sum(P(r)/r) over the distribution."""
    from afrc.polymer import PolymerObject
    po = PolymerObject(all_aa)
    dist, prob = po.get_end_to_end_distribution()
    numeric = np.sum(prob[1:] / dist[1:])
    assert po.get_mean_inverse_end_to_end_distance() == pytest.approx(numeric, rel=1e-3)


def test_mean_inverse_distance_is_not_inverse_of_mean(all_aa):
    """For a Gaussian chain <1/r> * <r> = 4/pi, not 1."""
    from afrc.polymer import PolymerObject
    po = PolymerObject(all_aa)
    product = po.get_mean_inverse_end_to_end_distance() * po.get_mean_end_to_end_distance('distribution')
    assert product == pytest.approx(4 / np.pi, rel=1e-3)


def test_mean_inverse_distance_zero_length_raises():
    from afrc.polymer import PolymerObject
    from afrc.exceptions import AFRCException
    with pytest.raises(AFRCException):
        PolymerObject('').get_mean_inverse_end_to_end_distance()


def test_prefactors_are_composition_weighted():
    """The scaling-law quantities use the composition-weighted per-residue prefactors."""
    from afrc.config import RIJ_R0, RIJ_RMS_R0, RG_R0

    seq = 'AAGGGK'
    po = PolymerObject(seq)
    n = len(seq)

    def weighted(table):
        return sum(table[aa] for aa in seq) / n

    assert po.get_mean_end_to_end_distance() == pytest.approx(weighted(RIJ_R0) * np.sqrt(n), rel=1e-12)
    assert po.get_mean_radius_of_gyration() == pytest.approx(weighted(RG_R0) * np.sqrt(n), rel=1e-12)
    assert po.RMS_Re_scaling == pytest.approx(weighted(RIJ_RMS_R0) * np.sqrt(n), rel=1e-12)


def test_end_to_end_distribution_has_the_scaling_law_rms(all_aa):
    po = PolymerObject(all_aa)
    dist, prob = po.get_end_to_end_distribution()
    assert np.sqrt(np.sum(prob * dist**2)) == pytest.approx(po.RMS_Re_scaling, rel=1e-4)


@pytest.mark.parametrize('n', [1, 20, 1000])
def test_distribution_grids_cover_the_tails(n):
    po = PolymerObject('G' * n)
    for dist, prob in (po.get_end_to_end_distribution(), po.get_radius_of_gyration_distribution()):
        assert prob[-1] < 1e-4 * prob.max()


def test_distributions_are_cached(all_aa):
    po = PolymerObject(all_aa)
    assert po.get_end_to_end_distribution()[1] is po.get_end_to_end_distribution()[1]
    assert po.get_radius_of_gyration_distribution()[1] is po.get_radius_of_gyration_distribution()[1]


def test_sampling_before_building_distributions(all_aa):
    """Sampling builds the distribution on demand."""
    np.random.seed(5)
    po = PolymerObject(all_aa)
    re = po.sample_end_to_end_distribution(dist_size=2000)
    rg = PolymerObject(all_aa).sample_radius_of_gyration_distribution(dist_size=2000)
    assert np.mean(re) == pytest.approx(po.get_mean_end_to_end_distance(), rel=0.05)
    assert np.mean(rg) == pytest.approx(po.get_mean_radius_of_gyration(), rel=0.05)
