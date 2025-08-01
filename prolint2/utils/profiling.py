"""
Performance benchmarking and profiling tools for ProLint2
"""

import time
import functools
import cProfile
import pstats
import io
from typing import Dict, Any, Optional, Callable, List, Tuple
import logging
import numpy as np
from contextlib import contextmanager

logger = logging.getLogger(__name__)


class PerformanceProfiler:
    """Profile and benchmark ProLint2 operations."""
    
    def __init__(self):
        self.timing_data = {}
        self.memory_data = {}
        self.operation_counts = {}
    
    def time_operation(self, operation_name: str):
        """Decorator to time operations."""
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                start_time = time.perf_counter()
                try:
                    result = func(*args, **kwargs)
                    return result
                finally:
                    end_time = time.perf_counter()
                    execution_time = end_time - start_time
                    
                    if operation_name not in self.timing_data:
                        self.timing_data[operation_name] = []
                        self.operation_counts[operation_name] = 0
                    
                    self.timing_data[operation_name].append(execution_time)
                    self.operation_counts[operation_name] += 1
                    
                    logger.debug(f"{operation_name} completed in {execution_time:.4f}s")
            
            return wrapper
        return decorator
    
    @contextmanager
    def profile_block(self, block_name: str):
        """Context manager to profile code blocks."""
        start_time = time.perf_counter()
        try:
            yield
        finally:
            end_time = time.perf_counter()
            execution_time = end_time - start_time
            
            if block_name not in self.timing_data:
                self.timing_data[block_name] = []
            
            self.timing_data[block_name].append(execution_time)
            logger.debug(f"Block '{block_name}' executed in {execution_time:.4f}s")
    
    def get_performance_summary(self) -> Dict[str, Any]:
        """Get summary of performance data."""
        summary = {}
        
        for operation, times in self.timing_data.items():
            times_array = np.array(times)
            summary[operation] = {
                'count': len(times),
                'total_time': times_array.sum(),
                'mean_time': times_array.mean(),
                'std_time': times_array.std(),
                'min_time': times_array.min(),
                'max_time': times_array.max(),
                'median_time': np.median(times_array)
            }
        
        return summary
    
    def print_performance_report(self):
        """Print formatted performance report."""
        summary = self.get_performance_summary()
        
        print("\n" + "="*80)
        print("PROLINT2 PERFORMANCE REPORT")
        print("="*80)
        print(f"{'Operation':<30} {'Count':<8} {'Total(s)':<10} {'Mean(s)':<10} {'Std(s)':<10}")
        print("-"*80)
        
        for operation, stats in summary.items():
            print(f"{operation:<30} {stats['count']:<8} "
                  f"{stats['total_time']:<10.4f} {stats['mean_time']:<10.4f} "
                  f"{stats['std_time']:<10.4f}")
        
        print("="*80)


def profile_function(func: Callable, *args, **kwargs) -> Tuple[Any, pstats.Stats]:
    """
    Profile a function call and return results with profiling stats.
    
    Parameters
    ----------
    func : Callable
        Function to profile
    *args, **kwargs
        Arguments to pass to the function
        
    Returns
    -------
    Tuple[Any, pstats.Stats]
        Function result and profiling statistics
    """
    profiler = cProfile.Profile()
    
    profiler.enable()
    try:
        result = func(*args, **kwargs)
    finally:
        profiler.disable()
    
    # Create stats object
    stats_stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stats_stream)
    
    return result, stats


def benchmark_contacts_computation(
    universe, 
    cutoffs: List[float] = [5.0, 7.0, 10.0, 15.0],
    n_runs: int = 3
) -> Dict[str, Dict[str, float]]:
    """
    Benchmark contact computation with different cutoffs.
    
    Parameters
    ----------
    universe : Universe
        ProLint2 Universe object
    cutoffs : List[float]
        List of cutoff distances to test
    n_runs : int
        Number of runs for each cutoff
        
    Returns
    -------
    Dict[str, Dict[str, float]]
        Benchmark results
    """
    results = {}
    
    for cutoff in cutoffs:
        times = []
        
        for run in range(n_runs):
            start_time = time.perf_counter()
            contacts = universe.compute_contacts(cutoff=cutoff)
            end_time = time.perf_counter()
            
            execution_time = end_time - start_time
            times.append(execution_time)
            
            # Clean up to ensure fair comparison
            del contacts
        
        times_array = np.array(times)
        results[f"cutoff_{cutoff}"] = {
            'mean_time': times_array.mean(),
            'std_time': times_array.std(),
            'min_time': times_array.min(),
            'max_time': times_array.max()
        }
        
        logger.info(f"Cutoff {cutoff}: {times_array.mean():.4f}±{times_array.std():.4f}s")
    
    return results


def estimate_computation_time(
    n_query_atoms: int,
    n_database_atoms: int, 
    n_frames: int,
    cutoff: float
) -> Dict[str, float]:
    """
    Estimate computation time based on system size.
    
    Parameters
    ----------
    n_query_atoms : int
        Number of query atoms
    n_database_atoms : int
        Number of database atoms
    n_frames : int
        Number of trajectory frames
    cutoff : float
        Distance cutoff
        
    Returns
    -------
    Dict[str, float]
        Time estimates
    """
    # Empirical scaling factors (these would need to be calibrated)
    base_time_per_frame = 1e-6  # seconds per atom pair per frame
    cutoff_scaling = (cutoff / 7.0) ** 2  # Quadratic scaling with cutoff
    
    # Rough estimates
    operations_per_frame = n_query_atoms * n_database_atoms
    estimated_time_per_frame = operations_per_frame * base_time_per_frame * cutoff_scaling
    total_estimated_time = estimated_time_per_frame * n_frames
    
    return {
        'estimated_time_per_frame': estimated_time_per_frame,
        'total_estimated_time': total_estimated_time,
        'operations_per_frame': operations_per_frame,
        'cutoff_scaling_factor': cutoff_scaling
    }


# Global profiler instance
global_profiler = PerformanceProfiler()


def benchmark_decorator(operation_name: str):
    """Decorator to automatically benchmark operations."""
    return global_profiler.time_operation(operation_name)
