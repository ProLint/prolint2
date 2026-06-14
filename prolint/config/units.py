"""Unit conversion and simulation parameters.

This module provides unit conversion factors and simulation parameter
containers for ProLint calculations.
"""

from enum import Enum
from dataclasses import dataclass, field


class UnitConversionFactor(Enum):
    """Conversion factors from time units to seconds.

    Each member's value represents the unit in seconds.

    Attributes
    ----------
    fs : float
        Femtoseconds (1e-15 s).
    ps : float
        Picoseconds (1e-12 s).
    ns : float
        Nanoseconds (1e-9 s).
    us : float
        Microseconds (1e-6 s).
    ms : float
        Milliseconds (1e-3 s).
    s : float
        Seconds (1.0 s).
    """

    fs = 1e-15
    ps = 1e-12
    ns = 1e-9
    us = 1e-6
    ms = 1e-3
    s = 1.0


@dataclass(frozen=True)
class SimulationParams:
    """Simulation parameters for contact calculations.

    Attributes
    ----------
    units : str
        Time units for output (default: "us").
    normalizer : str
        Normalization method (default: "actual_time").
    unit_conversion_factor : float
        Factor to convert trajectory time to output units.
    norm_factor : float
        Normalization factor for contact durations.
    """

    units: str = "us"
    normalizer: str = "actual_time"
    unit_conversion_factor: float = field(
        default_factory=lambda: UnitConversionFactor.ps.value
        / UnitConversionFactor.us.value
    )
    norm_factor: float = 1.0


# Legacy dict for backward compatibility
DEFAULT_SIM_PARAMS = {
    "units": "us",
    "normalizer": "actual_time",
    "unit_conversion_factor": UnitConversionFactor.ps.value
    / UnitConversionFactor.us.value,
    "norm_factor": 1,
}
