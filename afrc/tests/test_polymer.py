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
