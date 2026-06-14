"""Utility functions for ProLint.

This module provides optimized utility functions for contact
computation and data processing.
"""

import numpy as np


def fast_contiguous_segment_lengths(arr, multiplier: float = 1.0) -> np.ndarray:
    """Compute lengths of contiguous segments in a sorted array.

    Parameters
    ----------
    arr : array-like
        Sorted array of frame indices.
    multiplier : float, default=1.0
        Factor to multiply segment lengths by.

    Returns
    -------
    np.ndarray
        Array of segment lengths (contact durations).
    """
    if len(arr) == 0:
        return np.array([])

    # Calculate the differences between consecutive elements
    diffs = np.diff(arr)

    # Find the indices where the difference is greater than 1
    split_indices = np.where(diffs > 1)[0]

    # Calculate the segment lengths directly from the split_indices array using slicing
    segment_lengths = np.empty(split_indices.size + 1, dtype=int)
    if split_indices.size == 0:
        segment_lengths[0] = len(arr)
        return segment_lengths * multiplier
    segment_lengths[0] = split_indices[0] + 1
    segment_lengths[-1] = len(arr) - split_indices[-1] - 1
    segment_lengths[1:-1] = np.diff(split_indices)  # - 1

    return segment_lengths * multiplier


def fast_unique_comparison(residue_ids, database_ids, database_names):
    """Find unique residue-database pairs efficiently.

    Given parallel arrays of residue IDs, database IDs, and database names,
    returns the unique (residue_id, database_id) pairs with corresponding names.

    Parameters
    ----------
    residue_ids : np.ndarray
        Array of residue IDs.
    database_ids : np.ndarray
        Array of database molecule IDs.
    database_names : np.ndarray
        Array of database residue names.

    Returns
    -------
    tuple of np.ndarray
        (unique_residue_ids, unique_database_ids, unique_database_names)
    """
    # Handle empty input
    if len(residue_ids) == 0:
        return np.array([], dtype=residue_ids.dtype), np.array([], dtype=database_ids.dtype), np.array([], dtype=database_names.dtype)

    # Combine the arrays into a single 2D array
    combined_array = np.stack((residue_ids, database_ids), axis=-1)

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

    # Split the unique rows back into residue_ids and database_ids
    unique_residue_ids, unique_database_ids = unique_array[:, 0], unique_array[:, 1]

    # Extract the corresponding database_names using the sorted indices
    sorted_database_names = database_names[lex_sorted_indices]
    unique_database_names = sorted_database_names[unique_indices]

    return unique_residue_ids, unique_database_ids, unique_database_names
