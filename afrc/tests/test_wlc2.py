"""
Tests for ``afrc/polymer_models/wlc2.py`` - the O'Brien worm-like chain model.
"""

import numpy as np
import pytest

from afrc.polymer_models import wlc2


def test_distribution_is_valid_pmf(all_aa):
    model = wlc2.WormLikeChain2(all_aa)
    dist, prob = model.get_end_to_end_distribution()
    assert len(dist) == len(prob)
    assert np.all(prob >= 0)
    assert np.sum(prob) == pytest.approx(1.0)


def test_mean_rms_and_rg_positive(all_aa):
    model = wlc2.WormLikeChain2(all_aa)
    mean = model.get_mean_end_to_end_distance()
    rms = model.get_root_mean_squared_end_to_end_distance()
    rg = model.get_mean_radius_of_gyration()
    assert mean > 0
    assert rms >= mean
    assert rg > 0


def test_rg_matches_benoit_doty(all_aa):
    """Regression: the <Rg^2> expression previously had the wrong sign on the -lp^2 term."""
    lp, b = 3.0, 3.8
    N = len(all_aa)
    Lc = N * b

    model = wlc2.WormLikeChain2(all_aa, lp=lp, aa_size=b)

    # standard Benoit-Doty worm-like chain result
    rg_sq = Lc*lp/3 - lp**2 + 2*lp**3/Lc - 2*lp**4/Lc**2*(1 - np.exp(-Lc/lp))
    assert model.get_mean_radius_of_gyration() == pytest.approx(np.sqrt(rg_sq), rel=1e-9)


def test_rg_never_exceeds_the_rigid_rod_bound():
    """A chain can be no more extended than a rigid rod, so Rg^2 <= Lc^2/12.

    Regression: with the sign error on the -lp^2 term this bound was violated for
    short chains (e.g. a 2-residue chain gave Rg^2 = 21.1 against a rod bound of 4.8).
    """
    b = 3.8
    for n in (2, 3, 5, 10, 50, 200):
        seq = 'A' * n
        Lc = n * b
        rg = wlc2.WormLikeChain2(seq, lp=3.0, aa_size=b).get_mean_radius_of_gyration()
        assert rg**2 <= Lc**2/12


def test_rg_grows_with_persistence_length(all_aa):
    flexible = wlc2.WormLikeChain2(all_aa, lp=2.0).get_mean_radius_of_gyration()
    stiff = wlc2.WormLikeChain2(all_aa, lp=5.0).get_mean_radius_of_gyration()
    assert stiff > flexible


def test_exception_is_real_exception():
    assert issubclass(wlc2.WLC2Exception, Exception)


def test_rejects_chain_shorter_than_persistence_length():
    """The contour length (N * aa_size), not the residue count, is compared with lp."""
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2('AA', lp=20.0)

    # 2 residues is a contour length of 7.6 A, comfortably above the default lp = 3 A
    assert wlc2.WormLikeChain2('AA').get_mean_radius_of_gyration() > 0


def test_rejects_non_positive_lp(all_aa):
    with pytest.raises(wlc2.WLC2Exception):
        wlc2.WormLikeChain2(all_aa, lp=-1)
