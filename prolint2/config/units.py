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
