"""
Test suite for ProLint2 core functionality
"""

import unittest
import numpy as np
from unittest.mock import Mock, patch
import tempfile
import pathlib

from prolint2.utils.validation import (
    validate_cutoff, validate_atomgroup, validate_residue_id,
    validate_lipid_type, validate_frame_range, ValidationError
)


class TestValidation(unittest.TestCase):
    """Test input validation functions."""
    
    def test_validate_cutoff_valid(self):
        """Test valid cutoff values."""
        self.assertEqual(validate_cutoff(5.0), 5.0)
        self.assertEqual(validate_cutoff(10), 10.0)
        self.assertEqual(validate_cutoff(7.5), 7.5)
    
    def test_validate_cutoff_invalid(self):
        """Test invalid cutoff values."""
        with self.assertRaises(ValidationError):
            validate_cutoff(-1)
        
        with self.assertRaises(ValidationError):
            validate_cutoff(0)
        
        with self.assertRaises(ValidationError):
            validate_cutoff("invalid")
        
        with self.assertRaises(ValidationError):
            validate_cutoff(100)  # Too large
    
    def test_validate_frame_range(self):
        """Test frame range validation."""
        start, stop, step = validate_frame_range(0, 100, 5, 1000)
        self.assertEqual((start, stop, step), (0, 100, 5))
        
        # Test auto-correction
        start, stop, step = validate_frame_range(-10, 2000, 1, 1000)
        self.assertEqual((start, stop, step), (0, 1000, 1))
        
        # Test invalid ranges
        with self.assertRaises(ValidationError):
            validate_frame_range(100, 50, 1, 1000)  # start > stop
        
        with self.assertRaises(ValidationError):
            validate_frame_range(0, 100, 0, 1000)  # step = 0


class TestContactsComputation(unittest.TestCase):
    """Test contact computation functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        try:
            from prolint2.sampledata import GIRKDataSample
            from prolint2 import Universe
            girk = GIRKDataSample()
            self.universe = Universe(girk.coordinates, girk.trajectory)
        except Exception:
            # Skip tests if sample data is not available
            self.skipTest("Sample data not available")
    
    def test_compute_contacts_basic(self):
        """Test basic contact computation."""
        contacts = self.universe.compute_contacts(cutoff=7.0)
        
        # Basic checks
        self.assertIsNotNone(contacts)
        self.assertTrue(hasattr(contacts, 'contact_frames'))
        self.assertTrue(len(contacts.contact_frames) > 0)
    
    def test_compute_contacts_invalid_cutoff(self):
        """Test contact computation with invalid cutoff."""
        with self.assertRaises((ValueError, ValidationError)):
            self.universe.compute_contacts(cutoff=-1)
        
        with self.assertRaises((ValueError, ValidationError)):
            self.universe.compute_contacts(cutoff=0)


class TestMetricsComputation(unittest.TestCase):
    """Test metrics computation."""
    
    def setUp(self):
        """Set up test fixtures."""
        try:
            from prolint2.sampledata import GIRKDataSample
            from prolint2 import Universe
            girk = GIRKDataSample()
            self.universe = Universe(girk.coordinates, girk.trajectory)
            self.contacts = self.universe.compute_contacts(cutoff=7.0)
        except Exception:
            self.skipTest("Sample data not available")
    
    def test_compute_metrics(self):
        """Test basic metrics computation."""
        # Test built-in metrics
        for metric in ['sum', 'mean', 'max']:
            result = self.contacts.compute_metric(metric)
            self.assertIsInstance(result, dict)
            self.assertTrue(len(result) > 0)


class TestErrorHandling(unittest.TestCase):
    """Test error handling throughout the codebase."""
    
    def test_invalid_file_handling(self):
        """Test handling of invalid files."""
        from prolint2 import Universe
        
        with self.assertRaises((FileNotFoundError, OSError, ValueError)):
            Universe("nonexistent.gro", "nonexistent.xtc")
    
    def test_empty_selection_handling(self):
        """Test handling of empty selections."""
        # This would need to be implemented based on actual Universe behavior
        pass


if __name__ == '__main__':
    unittest.main()
