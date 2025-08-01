"""
General utility functions for ProLint2
======================================

This module contains general-purpose utility functions used throughout
the ProLint2 package.
"""

import numpy as np
from typing import Tuple


def fast_unique_comparison(
    residue_ids: np.ndarray, 
    lipid_ids: np.ndarray, 
    lipid_names: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get the unique combinations of residue and lipid ids using vectorized operations.

    This function efficiently finds unique combinations of residue-lipid pairs
    using lexicographic sorting, which is faster than naive approaches for
    large datasets.

    Parameters
    ----------
    residue_ids : np.ndarray
        Array of residue ids
    lipid_ids : np.ndarray
        Array of lipid ids  
    lipid_names : np.ndarray
        Array of lipid names corresponding to lipid_ids

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        Tuple containing:
        - unique_residue_ids: Array of unique residue ids
        - unique_lipid_ids: Array of unique lipid ids
        - unique_lipid_names: Array of corresponding unique lipid names

    Examples
    --------
    >>> residue_ids = np.array([1, 1, 2, 2, 1])
    >>> lipid_ids = np.array([10, 10, 20, 30, 10]) 
    >>> lipid_names = np.array(['POPC', 'POPC', 'POPE', 'CHOL', 'POPC'])
    >>> ur, ul, un = fast_unique_comparison(residue_ids, lipid_ids, lipid_names)
    >>> print(ur)  # [1 2 2]
    >>> print(ul)  # [10 20 30] 
    >>> print(un)  # ['POPC' 'POPE' 'CHOL']
    """
    # Input validation
    if len(residue_ids) != len(lipid_ids) != len(lipid_names):
        raise ValueError("All input arrays must have the same length")
    
    if len(residue_ids) == 0:
        return np.array([], dtype=residue_ids.dtype), \
               np.array([], dtype=lipid_ids.dtype), \
               np.array([], dtype=lipid_names.dtype)
    
    # Combine the arrays into a single 2D array for efficient sorting
    combined_array = np.stack((residue_ids, lipid_ids), axis=-1)

    # Get lexicographically sorted indices
    lex_sorted_indices = np.lexsort((combined_array[:, 1], combined_array[:, 0]))

    # Sort the combined array by the sorted indices
    sorted_array = combined_array[lex_sorted_indices]

    # Calculate row-wise differences between consecutive sorted rows
    row_diffs = np.diff(sorted_array, axis=0)

    # Find the indices where the differences are non-zero
    unique_indices = np.where(np.any(row_diffs != 0, axis=1))[0]

    # Add the first index (0) to unique_indices, as it's always unique
    unique_indices = np.concatenate(([0], unique_indices + 1))

    # Extract the unique rows using the indices
    unique_array = sorted_array[unique_indices]

    # Split the unique rows back into residue_ids and lipid_ids
    unique_residue_ids, unique_lipid_ids = unique_array[:, 0], unique_array[:, 1]

    # Extract the corresponding lipid_names using the sorted indices
    sorted_lipid_names = lipid_names[lex_sorted_indices]
    unique_lipid_names = sorted_lipid_names[unique_indices]

    return unique_residue_ids, unique_lipid_ids, unique_lipid_names