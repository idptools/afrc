"""
Shared pytest fixtures and constants for the afrc test suite.

Fixtures defined here are automatically available to every ``test_*.py`` module
in this directory without an explicit import.
"""

import pytest

import afrc


# a real, fully-disordered test sequence (FUS LC-like) used throughout
TEST_SEQ = (
    'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYG'
    'STGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQG'
)

# a short sequence containing all twenty amino acids exactly once
ALL_AA = 'ACDEFGHIKLMNPQRSTVWY'


@pytest.fixture(scope="session")
def test_seq():
    """The long disordered reference sequence used across the suite."""
    return TEST_SEQ


@pytest.fixture(scope="session")
def all_aa():
    """A short sequence containing each amino acid exactly once."""
    return ALL_AA


@pytest.fixture(scope="module")
def protein():
    """An ``AnalyticalFRC`` object built from the reference sequence."""
    return afrc.AnalyticalFRC(TEST_SEQ)
