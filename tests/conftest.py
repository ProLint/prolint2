"""Global pytest fixtures for ProLint testing.

This module provides reusable fixtures for testing ProLint components.
Fixtures are organized by category:
- MDAnalysis mock fixtures
- Contact data fixtures
- Analysis result fixtures
- Matplotlib mock fixtures
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch


# ============================================
# MDAnalysis Mock Fixtures
# ============================================


@pytest.fixture
def mock_residue():
    """Create a mock MDAnalysis Residue."""
    residue = MagicMock()
    residue.resid = 1
    residue.resname = "ALA"
    residue.segid = "PROA"
    return residue


@pytest.fixture
def mock_atom_group():
    """Create a mock MDAnalysis AtomGroup.

    Returns a MagicMock with properties matching a typical protein AtomGroup:
    - 100 atoms across 10 residues
    - Residue IDs 1-10
    - Mix of ALA and GLY residues
    """
    ag = MagicMock()
    ag.n_atoms = 100
    ag.n_residues = 10
    ag.residues.resids = np.arange(1, 11)
    ag.residues.resnames = np.array(["ALA"] * 5 + ["GLY"] * 5)
    ag.positions = np.random.rand(100, 3) * 50  # 50 Angstrom box
    ag.universe = MagicMock()
    return ag


@pytest.fixture
def mock_lipid_atom_group():
    """Create a mock MDAnalysis AtomGroup for lipids.

    Returns a MagicMock representing a membrane:
    - 500 atoms across 50 lipids
    - Mix of POPC, POPE, and CHOL
    """
    ag = MagicMock()
    ag.n_atoms = 500
    ag.n_residues = 50
    ag.residues.resids = np.arange(1, 51)
    ag.residues.resnames = np.array(["POPC"] * 20 + ["POPE"] * 20 + ["CHOL"] * 10)
    ag.positions = np.random.rand(500, 3) * 50
    ag.universe = MagicMock()
    return ag


@pytest.fixture
def mock_trajectory():
    """Create a mock MDAnalysis trajectory."""
    traj = MagicMock()
    traj.n_frames = 100
    traj.dt = 1.0  # 1 ps per frame
    traj.totaltime = 100.0  # 100 ps total
    traj.ts = MagicMock()
    traj.ts.frame = 0
    traj.ts.time = 0.0
    return traj


@pytest.fixture
def mock_mda_universe(mock_atom_group, mock_trajectory):
    """Create a mock MDAnalysis Universe.

    Combines mock atom group and trajectory into a complete universe mock.
    """
    u = MagicMock()
    u.trajectory = mock_trajectory
    u.atoms = mock_atom_group
    u.select_atoms = MagicMock(return_value=mock_atom_group)

    # Make the trajectory iterable
    def trajectory_iter():
        for i in range(mock_trajectory.n_frames):
            mock_trajectory.ts.frame = i
            mock_trajectory.ts.time = i * mock_trajectory.dt
            yield mock_trajectory.ts

    u.trajectory.__iter__ = trajectory_iter
    u.trajectory.__len__ = lambda: mock_trajectory.n_frames

    return u


@pytest.fixture
def mock_prolint_universe(mock_mda_universe, mock_atom_group, mock_lipid_atom_group):
    """Create a mock ProLint Universe with query and database groups.

    This fixture creates a MagicMock that behaves like a ProLint Universe
    without requiring actual file loading or MDAnalysis Universe initialization.
    """
    # Use a MagicMock to avoid MDAnalysis internal property setters
    u = MagicMock()
    u._query = mock_atom_group
    u._database = mock_lipid_atom_group
    u.params = {
        "units": "us",
        "normalizer": "counts",
        "norm_factor": 1.0,
        "unit_conversion_factor": 1.0,
    }
    u.trajectory = mock_mda_universe.trajectory
    u.atoms = mock_mda_universe.atoms
    u.select_atoms = mock_mda_universe.select_atoms

    # Set up query/database properties
    u.query = mock_atom_group
    u.database = mock_lipid_atom_group

    # Set up compute_contacts method
    def compute_contacts(**kwargs):
        mock_result = MagicMock()
        mock_result.contacts = {}
        mock_result.contact_frames = {}
        return mock_result

    u.compute_contacts = MagicMock(side_effect=compute_contacts)

    return u


# ============================================
# Contact Data Fixtures
# ============================================


@pytest.fixture
def sample_contact_frames():
    """Sample contact frame data structure.

    Format: {query_resid: {database_resid: [frame_indices]}}
    Represents which frames each query-database pair is in contact.
    """
    return {
        1: {10: [0, 1, 2, 5, 6, 7], 11: [3, 4, 8, 9]},
        2: {10: [0, 1, 2], 12: [5, 6, 7, 8, 9]},
        3: {11: list(range(100))},  # Full occupancy for testing
        4: {10: [0], 11: [1], 12: [2]},  # Minimal contacts
    }


@pytest.fixture
def sample_contact_durations():
    """Sample contact duration data structure.

    Format: {query_resid: {lipid_type: {database_resid: [durations]}}}
    Represents duration of each contact event in frames.
    """
    return {
        1: {"CHOL": {10: [3, 3], 11: [2, 2]}, "POPC": {20: [5]}},
        2: {"CHOL": {10: [3], 12: [5]}, "POPE": {30: [10, 5]}},
        3: {"POPC": {11: [100]}},  # Long-lived contact
        4: {"CHOL": {10: [1], 11: [1], 12: [1]}},  # Short contacts
    }


@pytest.fixture
def mock_computed_contacts(sample_contact_frames, sample_contact_durations):
    """Create a mock ComputedContacts object.

    Provides a minimal mock for testing analysis methods.
    """
    contacts = MagicMock()
    contacts.contact_frames = sample_contact_frames
    contacts.contacts = sample_contact_durations
    contacts.n_frames = 100
    contacts.query_resids = list(sample_contact_frames.keys())
    contacts.database_resnames = ["CHOL", "POPC", "POPE"]
    return contacts


# ============================================
# Analysis Result Fixtures
# ============================================


@pytest.fixture
def sample_timeseries_result():
    """Sample AnalysisResult for timeseries analysis."""
    from prolint.analysis.base import AnalysisResult

    n_residues = 10
    n_frames = 100

    return AnalysisResult(
        data={
            "query_residues": [
                {"resid": i, "resname": "ALA"} for i in range(1, n_residues + 1)
            ],
            "frames": list(range(n_frames)),
            "values": np.random.randint(0, 5, (n_residues, n_frames)).tolist(),
        },
        metadata={
            "analysis_type": "timeseries",
            "database_type": "CHOL",
            "n_frames": n_frames,
        },
    )


@pytest.fixture
def sample_metrics_result():
    """Sample AnalysisResult for metrics analysis."""
    from prolint.analysis.base import AnalysisResult

    n_residues = 10

    return AnalysisResult(
        data={
            "residues": [
                {"resid": i, "resname": "ALA"} for i in range(1, n_residues + 1)
            ],
            "values": np.random.rand(n_residues).tolist(),
            "metric": "occupancy",
        },
        metadata={
            "analysis_type": "metrics",
            "database_type": "CHOL",
            "metric": "occupancy",
        },
    )


@pytest.fixture
def sample_kinetics_result():
    """Sample AnalysisResult for kinetics analysis."""
    from prolint.analysis.base import AnalysisResult

    return AnalysisResult(
        data={
            "query_residue": {"resid": 42, "resname": "TRP"},
            "database_type": "CHOL",
            "kinetics": {
                "occupancy": 0.45,
                "n_events": 15,
                "mean_residence_time": 8.5,
                "koff": 0.118,
                "kon": 0.096,
            },
            "survival_times": list(range(1, 51)),
            "survival_curve": np.exp(-0.1 * np.arange(50)).tolist(),
            "residence_times": [3, 5, 8, 12, 5, 3, 10, 7, 9, 4, 6, 8, 11, 2, 7],
        },
        metadata={
            "analysis_type": "kinetics",
            "query_residue": 42,
            "database_type": "CHOL",
            "mode": "accumulated",
        },
    )


@pytest.fixture
def sample_density_result():
    """Sample AnalysisResult for density map analysis."""
    from prolint.analysis.base import AnalysisResult

    bins = 50

    return AnalysisResult(
        data={
            "density": np.random.rand(bins, bins).tolist(),
            "x_edges": np.linspace(-50, 50, bins + 1).tolist(),
            "y_edges": np.linspace(-50, 50, bins + 1).tolist(),
            "database_type": "CHOL",
        },
        metadata={
            "analysis_type": "density_map",
            "database_types": ["CHOL"],
            "bins": bins,
        },
    )


# ============================================
# Matplotlib Mock Fixtures
# ============================================


@pytest.fixture
def mock_matplotlib():
    """Mock matplotlib to avoid display issues in tests.

    Use this fixture for tests that call plotting functions.
    Prevents any GUI windows from opening during testing.
    """
    with patch("matplotlib.pyplot.figure") as mock_fig:
        fig_instance = MagicMock()
        ax_instance = MagicMock()
        fig_instance.add_subplot.return_value = ax_instance
        fig_instance.gca.return_value = ax_instance
        mock_fig.return_value = fig_instance
        yield mock_fig, fig_instance, ax_instance


@pytest.fixture(autouse=True)
def matplotlib_backend():
    """Set matplotlib to non-interactive backend for all tests.

    This runs automatically for all tests to prevent GUI issues.
    """
    import matplotlib

    matplotlib.use("Agg")


# ============================================
# Utility Fixtures
# ============================================


@pytest.fixture
def temp_pdb_file(tmp_path):
    """Create a minimal valid PDB file for testing.

    Returns path to a temporary PDB file with minimal valid content.
    """
    pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.246   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.986  -0.760   1.220  1.00  0.00           C
END
"""
    pdb_path = tmp_path / "test.pdb"
    pdb_path.write_text(pdb_content)
    return pdb_path


@pytest.fixture
def temp_gro_file(tmp_path):
    """Create a minimal valid GRO file for testing.

    Returns path to a temporary GRO file with minimal valid content.
    """
    gro_content = """\
Minimal GRO file
    5
    1ALA      N    1   0.000   0.000   0.000
    1ALA     CA    2   0.146   0.000   0.000
    1ALA      C    3   0.201   0.142   0.000
    1ALA      O    4   0.125   0.239   0.000
    1ALA     CB    5   0.199  -0.076   0.122
   5.00000   5.00000   5.00000
"""
    gro_path = tmp_path / "test.gro"
    gro_path.write_text(gro_content)
    return gro_path
