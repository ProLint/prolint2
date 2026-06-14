"""Unit tests for prolint.analysis module."""

import pytest
import numpy as np
from unittest.mock import MagicMock


class TestAnalysisResult:
    """Tests for AnalysisResult dataclass."""

    def test_creates_with_data_and_metadata(self):
        """AnalysisResult stores data and metadata."""
        from prolint.analysis.base import AnalysisResult

        data = {"values": [1, 2, 3]}
        metadata = {"type": "test"}
        result = AnalysisResult(data=data, metadata=metadata)

        assert result.data == data
        assert result.metadata == metadata


class TestBaseAnalysis:
    """Tests for BaseAnalysis abstract class."""

    def test_stores_universe_and_contacts(self):
        """BaseAnalysis stores universe and contacts references."""
        from prolint.analysis.base import BaseAnalysis

        class ConcreteAnalysis(BaseAnalysis):
            def run(self, **kwargs):
                return None

        mock_universe = MagicMock()
        mock_contacts = MagicMock()
        analysis = ConcreteAnalysis(mock_universe, mock_contacts)

        assert analysis.universe is mock_universe
        assert analysis.contacts is mock_contacts

    def test_rejects_none_inputs(self):
        """BaseAnalysis raises ValueError for None inputs."""
        from prolint.analysis.base import BaseAnalysis

        class ConcreteAnalysis(BaseAnalysis):
            def run(self, **kwargs):
                return None

        with pytest.raises(ValueError):
            ConcreteAnalysis(None, MagicMock())

        with pytest.raises(ValueError):
            ConcreteAnalysis(MagicMock(), None)

    def test_get_frame_range(self):
        """_get_frame_range returns correct frame ranges."""
        from prolint.analysis.base import BaseAnalysis

        class ConcreteAnalysis(BaseAnalysis):
            def run(self, **kwargs):
                return None

        mock_universe = MagicMock()
        mock_universe.trajectory.n_frames = 100
        analysis = ConcreteAnalysis(mock_universe, MagicMock())

        assert analysis._get_frame_range() == list(range(100))
        assert analysis._get_frame_range(frame_start=10, frame_end=50) == list(range(10, 50))
        assert analysis._get_frame_range(frame_start=0, frame_end=20, frame_step=5) == [0, 5, 10, 15]


class TestAnalysisRegistry:
    """Tests for AnalysisRegistry."""

    def test_all_nine_analyses_registered(self):
        """All 9 analysis types are registered."""
        from prolint.analysis import AnalysisRegistry

        available = AnalysisRegistry.available()
        expected = [
            "timeseries", "database_contacts", "kinetics", "density_map",
            "radial_density", "shared_contacts", "distances", "atom_distances", "metrics",
        ]

        assert len(available) == 9
        for name in expected:
            assert name in available

    def test_create_returns_instance(self):
        """create returns analysis instance with correct references."""
        from prolint.analysis import AnalysisRegistry

        mock_universe = MagicMock()
        mock_universe.trajectory.n_frames = 100
        mock_contacts = MagicMock()
        mock_contacts.contact_frames = {}

        analysis = AnalysisRegistry.create("timeseries", mock_universe, mock_contacts)

        assert analysis is not None
        assert analysis.universe is mock_universe

    def test_raises_for_unknown(self):
        """Registry raises ValueError for unknown analysis."""
        from prolint.analysis import AnalysisRegistry

        with pytest.raises(ValueError, match="Unknown analysis"):
            AnalysisRegistry.create("nonexistent", MagicMock(), MagicMock())

        with pytest.raises(ValueError, match="Unknown analysis"):
            AnalysisRegistry.get_class("nonexistent")
