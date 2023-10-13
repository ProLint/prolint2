from enum import Enum


# pylint: disable=invalid-name
class UnitConversionFactor(Enum):
    """Conversion factors for time units."""

    fs = 1e-15
    ps = 1e-12
    ns = 1e-9
    us = 1e-6
    ms = 1e-3
    s = 1.0


DEFAULT_SIM_PARAMS = {
    "units": "us",
    "normalizer": "actual_time",
    "unit_conversion_factor": UnitConversionFactor.ps.value
    / UnitConversionFactor.us.value,
    "norm_factor": 1,
}
