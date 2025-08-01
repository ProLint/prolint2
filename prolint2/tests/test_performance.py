"""
Performance benchmarks for ProLint2
"""

import time
from unittest.mock import Mock

import numpy as np
import pytest

from prolint2.utils.logging_config import setup_logger
from prolint2.utils.memory import MemoryMonitor, get_memory_usage
from prolint2.utils.validation import validate_cutoff, validate_frame_range


class TestValidationPerformance:
    """Performance tests for validation functions."""
    
    @pytest.mark.benchmark
    def test_validate_cutoff_performance(self, benchmark):
        """Benchmark cutoff validation."""
        result = benchmark(validate_cutoff, 7.0)
        assert result == 7.0
    
    @pytest.mark.benchmark
    def test_validate_frame_range_performance(self, benchmark):
        """Benchmark frame range validation."""
        result = benchmark(validate_frame_range, 0, 1000, 10, 10000)
        assert result == (0, 1000, 10)
    
    @pytest.mark.benchmark
    def test_multiple_validations_performance(self, benchmark):
        """Benchmark multiple validation calls."""
        def multiple_validations():
            for i in range(100):
                validate_cutoff(5.0 + i * 0.1)
        
        benchmark(multiple_validations)


class TestMemoryPerformance:
    """Performance tests for memory monitoring."""
    
    @pytest.mark.benchmark
    def test_memory_usage_performance(self, benchmark):
        """Benchmark memory usage retrieval."""
        result = benchmark(get_memory_usage)
        assert isinstance(result, float)
        assert result > 0
    
    @pytest.mark.benchmark
    def test_memory_monitor_performance(self, benchmark):
        """Benchmark memory monitor operations."""
        def monitor_operations():
            monitor = MemoryMonitor()
            monitor.start()
            for _ in range(10):
                monitor.update()
            return monitor.get_stats()
        
        result = benchmark(monitor_operations)
        assert isinstance(result, dict)


class TestLoggingPerformance:
    """Performance tests for logging system."""
    
    @pytest.mark.benchmark
    def test_logger_setup_performance(self, benchmark):
        """Benchmark logger setup."""
        def setup_logger_test():
            return setup_logger("perf_test", level="INFO")
        
        logger = benchmark(setup_logger_test)
        assert logger is not None
    
    @pytest.mark.benchmark
    def test_logging_messages_performance(self, benchmark):
        """Benchmark logging message performance."""
        logger = setup_logger("perf_test")
        
        def log_messages():
            for i in range(100):
                logger.info(f"Performance test message {i}")
        
        benchmark(log_messages)


class TestDataProcessingPerformance:
    """Performance tests for data processing operations."""
    
    @pytest.mark.benchmark
    def test_numpy_operations_performance(self, benchmark):
        """Benchmark numpy operations commonly used in ProLint2."""
        def numpy_operations():
            # Simulate distance calculations
            coords1 = np.random.random((1000, 3)) * 100
            coords2 = np.random.random((1000, 3)) * 100
            
            # Calculate distances
            distances = np.linalg.norm(coords1 - coords2, axis=1)
            
            # Find contacts within cutoff
            cutoff = 7.0
            contacts = distances < cutoff
            
            return np.sum(contacts)
        
        result = benchmark(numpy_operations)
        assert isinstance(result, (int, np.integer))
    
    @pytest.mark.benchmark
    def test_large_array_operations(self, benchmark):
        """Benchmark large array operations."""
        def large_array_ops():
            # Simulate large trajectory processing
            n_frames = 1000
            n_atoms = 5000
            
            # Create mock trajectory data
            trajectory = np.random.random((n_frames, n_atoms, 3)) * 100
            
            # Calculate center of mass for each frame
            com = np.mean(trajectory, axis=1)
            
            # Calculate distances from COM
            distances = np.linalg.norm(trajectory - com[:, np.newaxis, :], axis=2)
            
            return np.mean(distances)
        
        result = benchmark(large_array_ops)
        assert isinstance(result, (float, np.floating))


class TestMockPerformance:
    """Performance tests with mock objects."""
    
    @pytest.mark.benchmark
    def test_mock_universe_operations(self, benchmark):
        """Benchmark mock Universe operations."""
        def mock_universe_ops():
            # Create mock universe
            mock_universe = Mock()
            mock_universe.atoms = Mock()
            mock_universe.atoms.__len__ = Mock(return_value=10000)
            
            # Simulate operations
            n_atoms = len(mock_universe.atoms)
            mock_universe.select_atoms = Mock(return_value=Mock(__len__=Mock(return_value=n_atoms//2)))
            
            # Simulate contact calculation
            mock_contacts = Mock()
            mock_contacts.contact_frames = {i: {j: list(range(100)) for j in range(10)} for i in range(50)}
            mock_universe.compute_contacts = Mock(return_value=mock_contacts)
            
            contacts = mock_universe.compute_contacts(cutoff=7.0)
            return len(contacts.contact_frames)
        
        result = benchmark(mock_universe_ops)
        assert result == 50


if __name__ == '__main__':
    pytest.main([__file__, '--benchmark-only', '-v'])
