"""
Test configuration for ProLint2
"""

import pytest
import numpy as np
import tempfile
import pathlib
from unittest.mock import Mock, patch
import MDAnalysis as mda

from prolint2 import Universe
from prolint2.sampledata import GIRKDataSample


@pytest.fixture
def sample_universe():
    """Fixture providing a sample Universe for testing."""
    girk = GIRKDataSample()
    return Universe(girk.coordinates, girk.trajectory)


@pytest.fixture
def mock_trajectory():
    """Fixture providing a mock trajectory for testing."""
    with tempfile.NamedTemporaryFile(suffix='.xtc', delete=False) as f:
        yield f.name
    pathlib.Path(f.name).unlink(missing_ok=True)


@pytest.fixture
def mock_topology():
    """Fixture providing a mock topology for testing."""
    with tempfile.NamedTemporaryFile(suffix='.gro', delete=False) as f:
        yield f.name
    pathlib.Path(f.name).unlink(missing_ok=True)


class TestDataGenerator:
    """Helper class to generate test data."""
    
    @staticmethod
    def create_mock_atomgroup(n_atoms: int = 100) -> Mock:
        """Create a mock AtomGroup for testing."""
        mock_ag = Mock(spec=mda.AtomGroup)
        mock_ag.__len__ = Mock(return_value=n_atoms)
        mock_ag.n_atoms = n_atoms
        mock_ag.positions = np.random.random((n_atoms, 3)) * 100
        mock_ag.resids = np.arange(1, n_atoms + 1)
        mock_ag.resnames = ['LIG'] * n_atoms
        return mock_ag
    
    @staticmethod
    def create_test_contacts():
        """Create test contact data."""
        contacts = {}
        for resid in range(1, 11):
            contacts[resid] = {}
            for lipid_id in range(100, 110):
                contacts[resid][lipid_id] = list(np.random.randint(0, 1000, 50))
        return contacts
