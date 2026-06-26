"""
Tests for ``afrc/polymer.py`` - the internal ``PolymerObject`` class.
"""

import numpy as np
import pytest

from afrc.polymer import PolymerObject


def test_zero_length_polymer_object():
    po = PolymerObject('')
    assert po.zero_length is True
    assert np.all(po.sample_end_to_end_distribution(10) == 0.0)
    assert np.all(po.sample_radius_of_gyration_distribution(10) == 0.0)


def test_apparent_rms_bond_length_positive(all_aa):
    """Regression test: this method previously crashed (wrong attribute name)."""
    po = PolymerObject(all_aa)
    b = po.compute_apparent_rms_bond_length()
    assert b > 0
    assert np.isfinite(b)


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


def test_mean_radius_of_gyration_scaling_law_relationship(all_aa):
    po = PolymerObject(all_aa)
    re = po.get_mean_end_to_end_distance('scaling law')
    rg = po.get_mean_radius_of_gyration('scaling law')
    assert rg == pytest.approx(re / np.sqrt(6))


def test_mean_radius_of_gyration_distribution_positive(all_aa):
    po = PolymerObject(all_aa)
    assert po.get_mean_radius_of_gyration('distribution') > 0


def test_sampling_sizes(all_aa):
    po = PolymerObject(all_aa)
    assert len(po.sample_end_to_end_distribution(64)) == 64
    assert len(po.sample_radius_of_gyration_distribution(64)) == 64
