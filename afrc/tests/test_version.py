"""
Tests for version resolution in ``afrc/__init__.py``.
"""

import os
import sys
import types
from importlib import metadata

import pytest

import afrc


def test_version_is_a_non_empty_string():
    assert isinstance(afrc.__version__, str)
    assert len(afrc.__version__) > 0


def test_falls_back_to_package_metadata(monkeypatch):
    # a None entry in sys.modules makes `from ._version import ...` raise
    monkeypatch.setitem(sys.modules, 'afrc._version', None)
    monkeypatch.setattr(metadata, 'version', lambda name: '9.9.9')
    assert afrc._resolve_version() == '9.9.9'


def _raise_not_found(name):
    raise metadata.PackageNotFoundError(name)


def test_falls_back_to_versioningit_for_a_source_checkout(monkeypatch):
    seen = {}

    def fake_get_version(path):
        seen['path'] = path
        return '8.8.8'

    monkeypatch.setitem(sys.modules, 'afrc._version', None)
    monkeypatch.setattr(metadata, 'version', _raise_not_found)
    monkeypatch.setitem(sys.modules, 'versioningit', types.SimpleNamespace(get_version=fake_get_version))

    assert afrc._resolve_version() == '8.8.8'

    # versioningit is pointed at the repository root (the directory above the package)
    package_dir = os.path.dirname(afrc.__file__)
    assert os.path.normpath(seen['path']) == os.path.normpath(os.path.join(package_dir, os.pardir))


def test_fails_loudly_if_no_source_works(monkeypatch):
    def broken_get_version(path):
        raise RuntimeError('not a git repository')

    monkeypatch.setitem(sys.modules, 'afrc._version', None)
    monkeypatch.setattr(metadata, 'version', _raise_not_found)
    monkeypatch.setitem(sys.modules, 'versioningit', types.SimpleNamespace(get_version=broken_get_version))

    with pytest.raises(RuntimeError):
        afrc._resolve_version()
