"""
Input validation utilities for ProLint2
"""

import numpy as np
import MDAnalysis as mda
from typing import Union, List, Optional, Any
import os
import pathlib


class ValidationError(Exception):
    """Custom exception for validation errors"""
    pass


def validate_file_exists(filepath: Union[str, pathlib.Path]) -> pathlib.Path:
    """
    Validate that a file exists and is readable.
    
    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to the file
        
    Returns
    -------
    pathlib.Path
        Validated file path
        
    Raises
    ------
    ValidationError
        If file doesn't exist or is not readable
    """
    path = pathlib.Path(filepath)
    if not path.exists():
        raise ValidationError(f"File does not exist: {filepath}")
    if not path.is_file():
        raise ValidationError(f"Path is not a file: {filepath}")
    if not os.access(path, os.R_OK):
        raise ValidationError(f"File is not readable: {filepath}")
    return path


def validate_cutoff(cutoff: float) -> float:
    """
    Validate cutoff distance parameter.
    
    Parameters
    ----------
    cutoff : float
        Distance cutoff in Angstroms
        
    Returns
    -------
    float
        Validated cutoff
        
    Raises
    ------
    ValidationError
        If cutoff is invalid
    """
    if not isinstance(cutoff, (int, float)):
        raise ValidationError(f"Cutoff must be a number, got {type(cutoff)}")
    
    cutoff = float(cutoff)
    if cutoff <= 0:
        raise ValidationError(f"Cutoff must be positive, got {cutoff}")
    if cutoff > 50:  # Reasonable upper limit for MD simulations
        raise ValidationError(f"Cutoff seems unreasonably large: {cutoff} Å")
    
    return cutoff


def validate_atomgroup(atomgroup: Any, name: str = "AtomGroup") -> mda.AtomGroup:
    """
    Validate that an object is a valid MDAnalysis AtomGroup.
    
    Parameters
    ----------
    atomgroup : Any
        Object to validate
    name : str
        Name for error messages
        
    Returns
    -------
    mda.AtomGroup
        Validated AtomGroup
        
    Raises
    ------
    ValidationError
        If not a valid AtomGroup
    """
    if not isinstance(atomgroup, mda.AtomGroup):
        raise ValidationError(f"{name} must be an MDAnalysis AtomGroup, "
                            f"got {type(atomgroup)}")
    
    if len(atomgroup) == 0:
        raise ValidationError(f"{name} is empty")
    
    return atomgroup


def validate_residue_id(resid: int, atomgroup: mda.AtomGroup) -> int:
    """
    Validate that a residue ID exists in the given AtomGroup.
    
    Parameters
    ----------
    resid : int
        Residue ID to validate
    atomgroup : mda.AtomGroup
        AtomGroup to check against
        
    Returns
    -------
    int
        Validated residue ID
        
    Raises
    ------
    ValidationError
        If residue ID is invalid
    """
    if not isinstance(resid, int):
        raise ValidationError(f"Residue ID must be an integer, got {type(resid)}")
    
    available_resids = atomgroup.residues.resids
    if resid not in available_resids:
        raise ValidationError(f"Residue ID {resid} not found in AtomGroup. "
                            f"Available IDs: {available_resids[:10]}...")
    
    return resid


def validate_lipid_type(lipid_type: str, available_types: List[str]) -> str:
    """
    Validate lipid type name.
    
    Parameters
    ----------
    lipid_type : str
        Lipid type name
    available_types : List[str]
        List of available lipid types
        
    Returns
    -------
    str
        Validated lipid type
        
    Raises
    ------
    ValidationError
        If lipid type is invalid
    """
    if not isinstance(lipid_type, str):
        raise ValidationError(f"Lipid type must be a string, got {type(lipid_type)}")
    
    if lipid_type not in available_types:
        raise ValidationError(f"Lipid type '{lipid_type}' not found. "
                            f"Available types: {available_types}")
    
    return lipid_type


def validate_frame_range(start: int, stop: int, step: int, n_frames: int) -> tuple:
    """
    Validate frame range parameters.
    
    Parameters
    ----------
    start, stop, step : int
        Frame range parameters
    n_frames : int
        Total number of frames
        
    Returns
    -------
    tuple
        Validated (start, stop, step)
        
    Raises
    ------
    ValidationError
        If frame range is invalid
    """
    if not all(isinstance(x, int) for x in [start, stop, step]):
        raise ValidationError("Frame parameters must be integers")
    
    if step <= 0:
        raise ValidationError(f"Step must be positive, got {step}")
    
    if start < 0:
        start = 0
    
    if stop > n_frames:
        stop = n_frames
    
    if start >= stop:
        raise ValidationError(f"Invalid frame range: start={start}, stop={stop}")
    
    return start, stop, step
