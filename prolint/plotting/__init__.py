"""ProLint plotting module.

This module provides visualization functions and plotter classes for
contact analysis results. Use the :func:`plot` function for convenient
access to all plot types.

Examples
--------
>>> from prolint.plotting import plot
>>> fig, ax = plot("heatmap", result)
>>> fig, ax = plot("survival_curve", kinetics_result)
"""

# Base plotting infrastructure
from prolint.plotting.base import BasePlotter, PlottingRegistry, plot

# Theme
from prolint.plotting.theme import (
    COLORS,
    COLOR_SCALES,
    AMINO_ACID_COLORS,
    AMINO_ACID_ONE_LETTER,
    GRADIENTS,
    UNIT_LABELS,
    apply_prolint_style,
    get_colormap,
    get_color_for_value,
    interpolate_color,
    interpolate_gradient,
    get_unit_label,
    hex_to_rgb,
)

# Plotter classes
from prolint.plotting.heatmap import (
    HeatmapPlotter,
    DistanceHeatmapPlotter,
    DatabaseContactsHeatmapPlotter,
)
from prolint.plotting.network import NetworkPlotter, MAX_NETWORK_RESIDUES
from prolint.plotting.timeseries import (
    TimeSeriesPlotter,
    DistanceTimeSeriesPlotter,
    ContactEventsPlotter,
)
from prolint.plotting.kinetics import (
    SurvivalCurvePlotter,
    ResidenceDistributionPlotter,
)
from prolint.plotting.density import DensityMapPlotter, RadialDensityPlotter
from prolint.plotting.residues import ResidueMetricsPlotter, LogoGridPlotter

# Structure export
from prolint.plotting.structure import write_pdb

__all__ = [
    # Core
    "plot",
    "PlottingRegistry",
    "BasePlotter",
    # Plotters
    "HeatmapPlotter",
    "DistanceHeatmapPlotter",
    "DatabaseContactsHeatmapPlotter",
    "NetworkPlotter",
    "MAX_NETWORK_RESIDUES",
    "TimeSeriesPlotter",
    "DistanceTimeSeriesPlotter",
    "ContactEventsPlotter",
    "SurvivalCurvePlotter",
    "ResidenceDistributionPlotter",
    "DensityMapPlotter",
    "RadialDensityPlotter",
    "ResidueMetricsPlotter",
    "LogoGridPlotter",
    # Theme
    "COLORS",
    "COLOR_SCALES",
    "AMINO_ACID_COLORS",
    "AMINO_ACID_ONE_LETTER",
    "GRADIENTS",
    "UNIT_LABELS",
    "apply_prolint_style",
    "get_colormap",
    "get_color_for_value",
    "interpolate_color",
    "interpolate_gradient",
    "get_unit_label",
    "hex_to_rgb",
    # Structure
    "write_pdb",
]
