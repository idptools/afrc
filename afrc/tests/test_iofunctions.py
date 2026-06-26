"""
Tests for ``afrc/iofunctions.py`` - the ``validate_keyword`` helper.
"""

import pytest

from afrc.afrc import AFRCException
from afrc.iofunctions import validate_keyword


def test_valid_keyword_returned_lowercased():
    result = validate_keyword(['distribution', 'scaling law'], 'distribution', 'mode')
    assert result == 'distribution'


def test_keyword_is_case_insensitive():
    result = validate_keyword(['distribution', 'scaling law'], 'DISTRIBUTION', 'mode')
    assert result == 'distribution'


def test_invalid_keyword_raises_afrc_exception():
    with pytest.raises(AFRCException):
        validate_keyword(['distribution', 'scaling law'], 'nonsense', 'mode')


def test_non_string_keyword_raises_afrc_exception():
    """A non-string input cannot be lower-cased and should raise AFRCException."""
    with pytest.raises(AFRCException):
        validate_keyword(['distribution', 'scaling law'], 42, 'mode')
