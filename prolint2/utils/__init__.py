"""
ProLint2 Utilities
=================

This module provides utility functions for memory management, performance profiling,
input validation, and type definitions used throughout ProLint2.
"""

# Core utilities
from .utils import fast_unique_comparison

# Validation utilities
from .validation import (
    ValidationError,
    validate_file_exists,
    validate_cutoff,
    validate_atomgroup,
    validate_residue_id,
    validate_lipid_type,
    validate_frame_range
)

# Memory management
from .memory import (
    MemoryMonitor,
    memory_efficient_chunked_operation,
    monitor_memory,
    optimize_array_dtype,
    get_memory_usage
)

# Performance profiling
from .profiling import (
    PerformanceProfiler,
    profile_function,
    benchmark_contacts_computation,
    benchmark_decorator,
    global_profiler
)

# Logging configuration
from .logging_config import setup_logger, logger, get_logger

# Configuration management
from .config import Config, global_config, get_config

# Type definitions
from .types import (
    ResidueID, LipidID, AtomID, FrameIndex, ContactCount, MetricValue,
    Float64Array, Int64Array, BoolArray,
    ContactFrameDict, ResidueContactDict,
    LipidMetricDict, ResidueMetricDict, MetricResultDict,
    ConfigValue, ConfigDict,
    Position, Coordinates, TrajectoryFrame,
    AnalysisResult, PlotData
)

__all__ = [
    # Core utilities
    'fast_unique_comparison',
    
    # Validation
    'ValidationError', 'validate_file_exists', 'validate_cutoff',
    'validate_atomgroup', 'validate_residue_id', 'validate_lipid_type',
    'validate_frame_range',
    
    # Memory management
    'MemoryMonitor', 'memory_efficient_chunked_operation', 'monitor_memory',
    'optimize_array_dtype', 'get_memory_usage',
    
    # Performance profiling
    'PerformanceProfiler', 'profile_function', 'benchmark_contacts_computation',
    'benchmark_decorator', 'global_profiler',
    
    # Logging
    'setup_logger', 'logger', 'get_logger',
    
    # Configuration
    'Config', 'global_config', 'get_config'
    
    # Types
    'ResidueID', 'LipidID', 'AtomID', 'FrameIndex', 'ContactCount', 'MetricValue',
    'Float64Array', 'Int64Array', 'BoolArray',
    'ContactFrameDict', 'ResidueContactDict',
    'LipidMetricDict', 'ResidueMetricDict', 'MetricResultDict',
    'ConfigValue', 'ConfigDict',
    'Position', 'Coordinates', 'TrajectoryFrame',
    'AnalysisResult', 'PlotData'
]