"""
Test suite for ProLint2 validation utilities
"""

import pathlib
import tempfile
from unittest.mock import MagicMock, Mock

import MDAnalysis as mda
import numpy as np
import pytest

from prolint2.utils.validation import (
    ValidationError,
    validate_atomgroup,
    validate_cutoff,
    validate_file_exists,
    validate_frame_range,
    validate_lipid_type,
    validate_residue_id,
)


class TestFileValidation:
    """Test file validation functions."""
    
    def test_validate_file_exists_valid(self, tmp_path):
        """Test validation with existing file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")
        
        result = validate_file_exists(test_file)
        assert result == test_file
        assert isinstance(result, pathlib.Path)
    
    def test_validate_file_exists_nonexistent(self, tmp_path):
        """Test validation with non-existent file."""
        test_file = tmp_path / "nonexistent.txt"
        
        with pytest.raises(ValidationError, match="File does not exist"):
            validate_file_exists(test_file)
    
    def test_validate_file_exists_directory(self, tmp_path):
        """Test validation with directory instead of file."""
        test_dir = tmp_path / "test_dir"
        test_dir.mkdir()
        
        with pytest.raises(ValidationError, match="Path is not a file"):
            validate_file_exists(test_dir)
    
    def test_validate_file_exists_string_path(self, tmp_path):
        """Test validation with string path."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")
        
        result = validate_file_exists(str(test_file))
        assert result == test_file


class TestCutoffValidation:
    """Test cutoff validation."""
    
    @pytest.mark.parametrize("cutoff,expected", [
        (5.0, 5.0),
        (10, 10.0),
        (7.5, 7.5),
        (0.1, 0.1),
        (49.9, 49.9),
    ])
    def test_validate_cutoff_valid(self, cutoff, expected):
        """Test valid cutoff values."""
        assert validate_cutoff(cutoff) == expected
    
    @pytest.mark.parametrize("cutoff", [
        -1, 0, -5.5, 51, 100, "invalid", None, [], {}
    ])
    def test_validate_cutoff_invalid(self, cutoff):
        """Test invalid cutoff values."""
        with pytest.raises(ValidationError):
            validate_cutoff(cutoff)


class TestAtomGroupValidation:
    """Test AtomGroup validation."""
    
    def test_validate_atomgroup_valid(self):
        """Test with valid AtomGroup."""
        mock_ag = Mock(spec=mda.AtomGroup)
        mock_ag.__len__ = Mock(return_value=100)
        
        result = validate_atomgroup(mock_ag)
        assert result is mock_ag
    
    def test_validate_atomgroup_empty(self):
        """Test with empty AtomGroup."""
        mock_ag = Mock(spec=mda.AtomGroup)
        mock_ag.__len__ = Mock(return_value=0)
        
        with pytest.raises(ValidationError, match="is empty"):
            validate_atomgroup(mock_ag)
    
    def test_validate_atomgroup_invalid_type(self):
        """Test with invalid type."""
        with pytest.raises(ValidationError, match="must be an MDAnalysis AtomGroup"):
            validate_atomgroup("not an atomgroup")
    
    def test_validate_atomgroup_custom_name(self):
        """Test with custom name in error message."""
        with pytest.raises(ValidationError, match="CustomGroup must be"):
            validate_atomgroup("invalid", name="CustomGroup")


class TestResidueValidation:
    """Test residue ID validation."""
    
    def test_validate_residue_id_valid(self):
        """Test with valid residue ID."""
        mock_ag = Mock(spec=mda.AtomGroup)
        mock_residues = Mock()
        mock_residues.resids = np.array([1, 2, 3, 4, 5])
        mock_ag.residues = mock_residues
        
        result = validate_residue_id(3, mock_ag)
        assert result == 3
    
    def test_validate_residue_id_invalid(self):
        """Test with invalid residue ID."""
        mock_ag = Mock(spec=mda.AtomGroup)
        mock_residues = Mock()
        mock_residues.resids = np.array([1, 2, 3, 4, 5])
        mock_ag.residues = mock_residues
        
        with pytest.raises(ValidationError, match="Residue ID 10 not found"):
            validate_residue_id(10, mock_ag)
    
    def test_validate_residue_id_wrong_type(self):
        """Test with wrong type."""
        mock_ag = Mock(spec=mda.AtomGroup)
        
        with pytest.raises(ValidationError, match="must be an integer"):
            validate_residue_id("not_int", mock_ag)


class TestLipidTypeValidation:
    """Test lipid type validation."""
    
    def test_validate_lipid_type_valid(self):
        """Test with valid lipid type."""
        available = ["POPC", "POPE", "CHOL"]
        result = validate_lipid_type("POPC", available)
        assert result == "POPC"
    
    def test_validate_lipid_type_invalid(self):
        """Test with invalid lipid type."""
        available = ["POPC", "POPE", "CHOL"]
        
        with pytest.raises(ValidationError, match="not found"):
            validate_lipid_type("DPPC", available)
    
    def test_validate_lipid_type_wrong_type(self):
        """Test with wrong type."""
        available = ["POPC", "POPE", "CHOL"]
        
        with pytest.raises(ValidationError, match="must be a string"):
            validate_lipid_type(123, available)


class TestFrameRangeValidation:
    """Test frame range validation."""
    
    def test_validate_frame_range_valid(self):
        """Test with valid frame range."""
        start, stop, step = validate_frame_range(0, 100, 5, 1000)
        assert (start, stop, step) == (0, 100, 5)
    
    def test_validate_frame_range_auto_correction(self):
        """Test auto-correction of frame range."""
        # Negative start should be corrected to 0
        start, stop, step = validate_frame_range(-10, 100, 1, 1000)
        assert start == 0
        
        # Stop beyond n_frames should be corrected
        start, stop, step = validate_frame_range(0, 2000, 1, 1000)
        assert stop == 1000
    
    def test_validate_frame_range_invalid_step(self):
        """Test with invalid step."""
        with pytest.raises(ValidationError, match="Step must be positive"):
            validate_frame_range(0, 100, 0, 1000)
        
        with pytest.raises(ValidationError, match="Step must be positive"):
            validate_frame_range(0, 100, -1, 1000)
    
    def test_validate_frame_range_invalid_range(self):
        """Test with invalid range."""
        with pytest.raises(ValidationError, match="Invalid frame range"):
            validate_frame_range(100, 50, 1, 1000)
    
    def test_validate_frame_range_wrong_types(self):
        """Test with wrong types."""
        with pytest.raises(ValidationError, match="must be integers"):
            validate_frame_range("0", 100, 1, 1000)
        
        with pytest.raises(ValidationError, match="must be integers"):
            validate_frame_range(0, 100.5, 1, 1000)


class TestValidationIntegration:
    """Integration tests for validation functions."""
    
    def test_validation_error_hierarchy(self):
        """Test that ValidationError is properly defined."""
        assert issubclass(ValidationError, Exception)
        
        error = ValidationError("test message")
        assert str(error) == "test message"
    
    def test_validation_chain(self, tmp_path):
        """Test chaining multiple validations."""
        # Create a real file for testing
        test_file = tmp_path / "test.gro"
        test_file.write_text("test content")
        
        # Chain validations
        validated_file = validate_file_exists(test_file)
        validated_cutoff = validate_cutoff(7.0)
        
        assert validated_file.exists()
        assert validated_cutoff == 7.0
