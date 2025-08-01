"""
Type definitions for ProLint2
"""

from typing import Any, Dict, Generic, List, Optional, Tuple, TypeVar, Union

import numpy as np
from numpy.typing import NDArray

# Type aliases for better code readability
ResidueID = int
LipidID = int
AtomID = int
FrameIndex = int
ContactCount = int
MetricValue = float

# Array types
Float64Array = NDArray[np.float64]
Int64Array = NDArray[np.int64]
BoolArray = NDArray[np.bool_]

# Contact data structures
ContactFrameDict = Dict[LipidID, List[FrameIndex]]
ResidueContactDict = Dict[ResidueID, ContactFrameDict]

# Metric data structures
LipidMetricDict = Dict[LipidID, MetricValue]
ResidueMetricDict = Dict[str, LipidMetricDict]  # str is lipid type name
MetricResultDict = Dict[ResidueID, ResidueMetricDict]

# Configuration types
ConfigValue = Union[str, int, float, bool, List[str]]
ConfigDict = Dict[str, ConfigValue]

# Coordinate types
Position = Tuple[float, float, float]
Coordinates = NDArray[np.float64]  # Shape: (n_atoms, 3)
TrajectoryFrame = NDArray[np.float64]  # Shape: (n_atoms, 3)

# Analysis types
AnalysisResult = Dict[str, Any]
PlotData = Dict[str, Union[List, NDArray]]

# Generic types for extensibility
T = TypeVar('T')
AnalysisData = Generic[T]
