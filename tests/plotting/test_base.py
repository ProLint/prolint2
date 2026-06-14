"""Unit tests for prolint.plotting module."""

import pytest
import matplotlib
matplotlib.use('Agg')

import numpy as np
from unittest.mock import MagicMock, patch


class TestPlottingRegistry:
    """Tests for PlottingRegistry."""

    def test_all_thirteen_plotters_registered(self):
        """All 13 plotter types are registered."""
        from prolint.plotting import PlottingRegistry

        available = PlottingRegistry.available()
        expected = [
            "heatmap", "distance_heatmap", "database_contacts_heatmap",
            "timeseries", "distance_timeseries", "contact_events",
            "survival_curve", "residence_distribution",
            "density_map", "radial_density", "network",
            "residue_metrics", "logo_grid",
        ]

        assert len(available) == 13
        for name in expected:
            assert name in available

    def test_get_class_returns_plotter(self):
        """get_class returns plotter class."""
        from prolint.plotting import PlottingRegistry
        from prolint.plotting.base import BasePlotter

        cls = PlottingRegistry.get_class("heatmap")
        assert issubclass(cls, BasePlotter)

    def test_raises_for_unknown(self):
        """Registry raises ValueError for unknown plotter."""
        from prolint.plotting import PlottingRegistry

        with pytest.raises(ValueError):
            PlottingRegistry.get_class("nonexistent")

        with pytest.raises(ValueError):
            PlottingRegistry.plot("nonexistent", MagicMock())


class TestPlotFunction:
    """Tests for plot convenience function."""

    def test_plot_creates_figure(self):
        """plot function creates figure and axes."""
        from prolint.plotting import plot
        from prolint.analysis.base import AnalysisResult
        import matplotlib.pyplot as plt

        result = AnalysisResult(
            data={
                "frames": list(range(10)),
                "contact_counts": {1: list(range(10))},
            }
        )

        fig, ax = plot("heatmap", result)

        assert fig is not None
        assert ax is not None
        plt.close(fig)


class TestPlotterValidation:
    """Tests for plotter validation."""

    def test_heatmap_validates_timeseries_format(self):
        """HeatmapPlotter accepts timeseries format."""
        from prolint.plotting.heatmap import HeatmapPlotter
        from prolint.analysis.base import AnalysisResult

        valid = AnalysisResult(data={"frames": [0, 1], "contact_counts": {1: [1, 2]}})
        HeatmapPlotter.validate_result(valid)  # Should not raise

        invalid = AnalysisResult(data={"invalid": "data"})
        with pytest.raises(ValueError):
            HeatmapPlotter.validate_result(invalid)

    def test_timeseries_validates_required_keys(self):
        """TimeSeriesPlotter validates required keys."""
        from prolint.plotting.timeseries import TimeSeriesPlotter
        from prolint.analysis.base import AnalysisResult

        valid = AnalysisResult(data={"frames": [0, 1], "contact_counts": {1: [1, 2]}})
        TimeSeriesPlotter.validate_result(valid)  # Should not raise

        invalid = AnalysisResult(data={"frames": [0, 1]})  # Missing contact_counts
        with pytest.raises(ValueError):
            TimeSeriesPlotter.validate_result(invalid)

    def test_network_validates_matrix_and_labels(self):
        """NetworkPlotter validates matrix and labels."""
        from prolint.plotting.network import NetworkPlotter
        from prolint.analysis.base import AnalysisResult

        valid = AnalysisResult(data={"matrix": [[0, 1], [1, 0]], "labels": [1, 2]})
        NetworkPlotter.validate_result(valid)  # Should not raise

        invalid = AnalysisResult(data={"matrix": [[0, 1], [1, 0]]})  # Missing labels
        with pytest.raises(ValueError):
            NetworkPlotter.validate_result(invalid)
