"""ProLint configuration module.

This module provides configuration classes, enums, and utilities
for customizing ProLint behavior.
"""

from prolint.config.enums import TimeUnit, NormalizationMethod
from prolint.config.units import UnitConversionFactor, SimulationParams
from prolint.config.colors import (
    load_theme,
    COLORS,
    COLOR_SCALES,
    GRADIENTS,
    AMINO_ACID_COLORS,
    AMINO_ACID_ONE_LETTER,
    UNIT_LABELS,
    hex_to_rgb,
    color_to_tuple,
    interpolate_gradient,
    get_color_for_value,
    get_unit_label,
)
from prolint.config.logging import setup_logging, get_logger

__all__ = [
    "TimeUnit",
    "NormalizationMethod",
    "UnitConversionFactor",
    "SimulationParams",
    "load_theme",
    "COLORS",
    "COLOR_SCALES",
    "GRADIENTS",
    "AMINO_ACID_COLORS",
    "AMINO_ACID_ONE_LETTER",
    "UNIT_LABELS",
    "hex_to_rgb",
    "color_to_tuple",
    "interpolate_gradient",
    "get_color_for_value",
    "get_unit_label",
    "setup_logging",
    "get_logger",
]
