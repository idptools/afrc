"""
Unit and regression test for the afrc package.
"""

# Import package, test suite, and other packages as needed
import afrc
import pytest
import sys

def test_afrc_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "afrc" in sys.modules


def test_mean_rg():
    s = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQG'
    
    protein = afrc.AnalyticalFRC(s)

    # check mean rg code
    assert protein.get_mean_radius_of_gyration() - 30.788294820233368 < 0.001
    assert protein.get_mean_radius_of_gyration(calculation_mode='scaling law') - 29.341888777603817 < 0.001
    assert protein.get_mean_radius_of_gyration(calculation_mode='distribution') - 30.788294820233368 < 0.001
