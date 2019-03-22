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
