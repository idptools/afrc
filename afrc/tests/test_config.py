"""
Tests for ``afrc/config.py`` - the parameter lookup tables.
"""

from afrc.config import (
    AA_list,
    RIJ_RMS_R0,
    RIJ_R0,
    RG_X0,
    RG_R0,
    P_OF_R_RESOLUTION,
)


def test_aa_list_has_twenty_unique_residues():
    assert len(AA_list) == 20
    assert len(set(AA_list)) == 20


def test_resolution_is_positive():
    assert P_OF_R_RESOLUTION > 0


def test_lookup_tables_cover_all_amino_acids():
    for table in (RIJ_RMS_R0, RIJ_R0, RG_X0, RG_R0):
        assert set(table.keys()) == set(AA_list)


def test_lookup_table_values_are_positive():
    for table in (RIJ_RMS_R0, RIJ_R0, RG_X0, RG_R0):
        assert all(v > 0 for v in table.values())


def test_rg_r0_is_inverse_relationship_of_x0():
    # documented relationship in config.py: RG_R0 = 1.3858 * (1 / RG_X0)
    for aa in AA_list:
        assert abs(RG_R0[aa] - 1.3858 / RG_X0[aa]) < 1e-2
