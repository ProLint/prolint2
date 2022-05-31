"""
Unit and regression test for the ufcc package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pandas as pd
from pandas.util.testing import assert_frame_equal

import ufcc


def test_ufcc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ufcc" in sys.modules


def test_ufcc_initialize():
    """Initialize test, will always pass so long as the UFCC object can be initialized."""
    assert ufcc.UFCC('ufcc/data/coordinates.gro', 'ufcc/data/trajectory.xtc')


def test_contacts():
    """Contacts test, will always pass so long as the UFCC object can be initialized."""
    target_system = ufcc.UFCC('ufcc/data/coordinates.gro', 'ufcc/data/trajectory.xtc')
    target_system.contacts.compute()
    target_system.contacts.count_contacts()
    ref = pd.read_pickle('ufcc/data/counted_contacts.pkl')
    assert_frame_equal(target_system.contacts.counts, ref) 


