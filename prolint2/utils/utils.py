import numpy as np


def fast_unique_comparison(residue_ids, lipid_ids, lipid_names):
    """
    Get the unique combinations of residue and lipid ids. Vectorized implementation.

    Parameters
    ----------
    residue_ids : np.ndarray
        Array of residue ids.
    lipid_ids : np.ndarray
        Array of lipid ids.
    lipid_names : np.ndarray
        Array of lipid names.

    Returns
    -------
    np.ndarray: Array of unique combinations of residue and lipid ids.
    """
    # Combine the arrays into a single 2D array
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
