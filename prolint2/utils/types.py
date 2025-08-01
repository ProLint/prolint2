"""
Type definitions for ProLint2
=============================

This module provides type aliases and type definitions used throughout
ProLint2 for better code readability and type safety.
"""

from typing import Any, Dict, Generic, List, Optional, Tuple, TypeVar, Union
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

# Basic type aliases for better code readability
ResidueID = int
LipidID = int
AtomID = int
FrameIndex = int
ContactCount = int
MetricValue = float
FilePath = Union[str, Path]

# Array types - commonly used numpy array types
Float64Array = NDArray[np.float64]
Float32Array = NDArray[np.float32]
Int64Array = NDArray[np.int64]
Int32Array = NDArray[np.int32]
BoolArray = NDArray[np.bool_]
StringArray = NDArray[np.str_]

# Contact data structures
ContactFrameDict = Dict[LipidID, List[FrameIndex]]
ResidueContactDict = Dict[ResidueID, ContactFrameDict]

# Metric data structures  
LipidMetricDict = Dict[LipidID, MetricValue]
ResidueMetricDict = Dict[str, LipidMetricDict]  # str is lipid type name
MetricResultDict = Dict[ResidueID, ResidueMetricDict]

# Configuration types
ConfigValue = Union[str, int, float, bool, List[str], Dict[str, Any]]
ConfigDict = Dict[str, ConfigValue]

# Coordinate types
Position = Tuple[float, float, float]
Coordinates = NDArray[np.float64]  # Shape: (n_atoms, 3)
TrajectoryFrame = NDArray[np.float64]  # Shape: (n_atoms, 3)

# Analysis types
AnalysisResult = Dict[str, Any]
PlotData = Dict[str, Union[List, NDArray, float, int, str]]

# Time series data
TimeSeriesData = NDArray[np.float64]  # 1D array of values over time
TimePoint = float

# Distance and cutoff types
Distance = float
Cutoff = float
DistanceMatrix = NDArray[np.float64]  # Shape: (n_atoms, n_atoms)

# Generic types for extensibility
T = TypeVar('T')
AnalysisData = Generic[T]

# File format types
TrajectoryFile = FilePath
TopologyFile = FilePath
OutputFile = FilePath

# Selection strings (MDAnalysis style)
SelectionString = str

# Callback function types
ProgressCallback = Optional[callable]
ValidationCallback = Optional[callable]
