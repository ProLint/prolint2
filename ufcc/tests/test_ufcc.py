"""
Unit and regression test for the ufcc package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import ufcc


def test_ufcc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ufcc" in sys.modules

def test_ufcc_initialize():
    """Initialize test, will always pass so long as the UFCC object can be initialized."""
    assert ufcc.UFCC('ufcc/data/coordinates.gro', 'ufcc/data/trajectory.xtc')
