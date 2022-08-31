"""
Unit tests for the ufcc package.
"""

# Import package, test suite, and other packages as needed
import sys
import os
import pytest

import pandas as pd
from pandas.util.testing import assert_frame_equal

import ufcc
from ufcc.sampledata import GIRK


def test_ufcc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ufcc" in sys.modules


def test_ufcc_initialize():
    """Initialize test, will always pass so long as the UFCC object can be initialized."""
    assert ufcc.UFCC(GIRK.coordinates, GIRK.trajectory)

def test_serial_backend():
    """Serial backend test"""
    target_system = ufcc.UFCC(GIRK.coordinates, GIRK.trajectory)
    target_system.contacts.runner.backend = 'serial'
    target_system.contacts.compute()
    target_system.contacts.count_contacts()
    ref = pd.read_pickle(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'counted_contacts.pkl'))
    assert_frame_equal(target_system.contacts.counts, ref)

# def test_parallel_backend():
#     """Parallel backend test"""
#     target_system = ufcc.UFCC(GIRK.coordinates, GIRK.trajectory)
#     target_system.contacts.runner.backend = 'parallel'
#     target_system.contacts.runner.n_jobs = -1
#     target_system.contacts.compute()
#     target_system.contacts.count_contacts()
#     ref = pd.read_pickle('counted_contacts.pkl')
#     assert_frame_equal(target_system.contacts.counts, ref)

def test_ufcc_get_metrics():
    """Get metrics test"""
    target_system = ufcc.UFCC(GIRK.coordinates, GIRK.trajectory)
    target_system.contacts.runner.backend = 'serial'
    target_system.contacts.compute()
    target_system.contacts.count_contacts()
    target_system.contacts.get_metrics() 
    assert isinstance(target_system.contacts.contact_metrics, pd.DataFrame)