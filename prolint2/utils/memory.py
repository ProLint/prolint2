"""
Memory management utilities for ProLint2
"""

import psutil
import gc
import numpy as np
from typing import Optional, Dict, Any
import logging
from functools import wraps
import warnings

logger = logging.getLogger(__name__)


class MemoryMonitor:
    """Monitor and manage memory usage during computations."""
    
    def __init__(self, warning_threshold: float = 0.8, error_threshold: float = 0.95):
        """
        Initialize memory monitor.
        
        Parameters
        ----------
        warning_threshold : float
            Memory usage fraction to trigger warnings (default: 0.8)
        error_threshold : float
            Memory usage fraction to trigger errors (default: 0.95)
        """
        self.warning_threshold = warning_threshold
        self.error_threshold = error_threshold
        self.initial_memory = self.get_memory_usage()
    
    def get_memory_usage(self) -> Dict[str, float]:
        """Get current memory usage statistics."""
        process = psutil.Process()
        memory_info = process.memory_info()
        virtual_memory = psutil.virtual_memory()
        
        return {
            'rss_mb': memory_info.rss / 1024 / 1024,  # Resident set size
            'vms_mb': memory_info.vms / 1024 / 1024,  # Virtual memory size
            'percent': process.memory_percent(),
            'available_mb': virtual_memory.available / 1024 / 1024,
            'total_mb': virtual_memory.total / 1024 / 1024,
            'system_percent': virtual_memory.percent / 100.0
        }
    
    def check_memory(self, operation_name: str = "operation") -> None:
        """
        Check current memory usage and warn/error if thresholds exceeded.
        
        Parameters
        ----------
        operation_name : str
            Name of the operation being performed
            
        Raises
        ------
        MemoryError
            If memory usage exceeds error threshold
        """
        memory_stats = self.get_memory_usage()
        usage_fraction = memory_stats['system_percent']
        
        if usage_fraction >= self.error_threshold:
            raise MemoryError(
                f"Memory usage too high during {operation_name}: "
                f"{usage_fraction:.1%} (threshold: {self.error_threshold:.1%})"
            )
        elif usage_fraction >= self.warning_threshold:
            warnings.warn(
                f"High memory usage during {operation_name}: "
                f"{usage_fraction:.1%} (threshold: {self.warning_threshold:.1%})",
                UserWarning
            )
            logger.warning(f"High memory usage: {memory_stats}")
    
    def suggest_optimization(self, array_size_mb: float) -> Dict[str, Any]:
        """
        Suggest optimizations based on array size and available memory.
        
        Parameters
        ----------
        array_size_mb : float
            Size of array in MB
            
        Returns
        -------
        Dict[str, Any]
            Optimization suggestions
        """
        memory_stats = self.get_memory_usage()
        available_mb = memory_stats['available_mb']
        
        suggestions = {
            'use_chunking': array_size_mb > available_mb * 0.5,
            'suggested_chunk_size': max(1, int(available_mb * 0.2)),
            'use_memory_mapping': array_size_mb > 1000,  # > 1GB
            'reduce_precision': array_size_mb > available_mb * 0.3,
            'available_memory_mb': available_mb,
            'array_size_mb': array_size_mb
        }
        
        return suggestions


def memory_efficient_chunked_operation(
    data: np.ndarray, 
    operation_func: callable,
    chunk_size: Optional[int] = None,
    axis: int = 0
) -> np.ndarray:
    """
    Perform operations on large arrays in chunks to manage memory.
    
    Parameters
    ----------
    data : np.ndarray
        Input data array
    operation_func : callable
        Function to apply to each chunk
    chunk_size : int, optional
        Size of chunks. If None, automatically determined.
    axis : int
        Axis along which to chunk
        
    Returns
    -------
    np.ndarray
        Result of operation
    """
    monitor = MemoryMonitor()
    
    if chunk_size is None:
        # Estimate reasonable chunk size based on available memory
        data_size_mb = data.nbytes / 1024 / 1024
        suggestions = monitor.suggest_optimization(data_size_mb)
        
        if suggestions['use_chunking']:
            chunk_size = min(
                data.shape[axis] // 4,  # At least 4 chunks
                suggestions['suggested_chunk_size'] * 1024 * 1024 // data.itemsize
            )
            chunk_size = max(1, chunk_size)
        else:
            return operation_func(data)
    
    # Process in chunks
    results = []
    n_chunks = int(np.ceil(data.shape[axis] / chunk_size))
    
    logger.info(f"Processing {data.shape} array in {n_chunks} chunks of size {chunk_size}")
    
    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, data.shape[axis])
        
        # Create slice object for the chunk
        chunk_slice = [slice(None)] * data.ndim
        chunk_slice[axis] = slice(start_idx, end_idx)
        chunk = data[tuple(chunk_slice)]
        
        # Monitor memory before processing chunk
        monitor.check_memory(f"chunk {i+1}/{n_chunks}")
        
        # Process chunk
        result = operation_func(chunk)
        results.append(result)
        
        # Force garbage collection between chunks
        del chunk
        gc.collect()
    
    # Combine results
    return np.concatenate(results, axis=axis)


def monitor_memory(func):
    """Decorator to monitor memory usage of functions."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        monitor = MemoryMonitor()
        initial_memory = monitor.get_memory_usage()
        
        logger.debug(f"Starting {func.__name__}: {initial_memory}")
        
        try:
            result = func(*args, **kwargs)
            return result
        finally:
            final_memory = monitor.get_memory_usage()
            memory_diff = final_memory['rss_mb'] - initial_memory['rss_mb']
            
            logger.debug(
                f"Finished {func.__name__}: "
                f"Memory change: {memory_diff:+.1f} MB, "
                f"Current usage: {final_memory['rss_mb']:.1f} MB"
            )
    
    return wrapper


def optimize_array_dtype(array: np.ndarray, preserve_precision: bool = True) -> np.ndarray:
    """
    Optimize array dtype to reduce memory usage.
    
    Parameters
    ----------
    array : np.ndarray
        Input array
    preserve_precision : bool
        Whether to preserve numerical precision
        
    Returns
    -------
    np.ndarray
        Array with optimized dtype
    """
    if array.dtype == np.float64 and not preserve_precision:
        # Check if we can safely downcast to float32
        if np.allclose(array, array.astype(np.float32), rtol=1e-6):
            logger.info("Downcasting float64 to float32 to save memory")
            return array.astype(np.float32)
    
    elif array.dtype in [np.int64]:
        # Check if we can use smaller integer types
        min_val, max_val = array.min(), array.max()
        
        if np.iinfo(np.int32).min <= min_val <= max_val <= np.iinfo(np.int32).max:
            logger.info("Downcasting int64 to int32 to save memory")
            return array.astype(np.int32)
        elif np.iinfo(np.int16).min <= min_val <= max_val <= np.iinfo(np.int16).max:
            logger.info("Downcasting int64 to int16 to save memory")
            return array.astype(np.int16)
    
    return array
