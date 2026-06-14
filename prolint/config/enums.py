"""ProLint enumeration types.

This module provides enum classes for time units and normalization methods.
"""

from enum import Enum


class TimeUnit(str, Enum):
    """Time unit options for contact duration analysis.

    Attributes
    ----------
    FEMTOSECOND : str
        Femtoseconds ("fs").
    PICOSECOND : str
        Picoseconds ("ps").
    NANOSECOND : str
        Nanoseconds ("ns").
    MICROSECOND : str
        Microseconds ("us").
    MILLISECOND : str
        Milliseconds ("ms").
    SECOND : str
        Seconds ("s").
    """

    FEMTOSECOND = "fs"
    PICOSECOND = "ps"
    NANOSECOND = "ns"
    MICROSECOND = "us"
    MILLISECOND = "ms"
    SECOND = "s"

    def __str__(self) -> str:
        """Return the string value of the enum."""
        return self.value


class NormalizationMethod(str, Enum):
    """Normalization methods for contact durations.

    Attributes
    ----------
    COUNTS : str
        Report durations as frame counts ("counts").
    ACTUAL_TIME : str
        Report durations in time units ("actual_time").
    """

    COUNTS = "counts"
    ACTUAL_TIME = "actual_time"

    def __str__(self) -> str:
        """Return the string value of the enum."""
        return self.value
