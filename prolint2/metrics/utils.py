from typing import Iterable, List, Dict
from itertools import chain

import numpy as np


def fast_filter_resids_by_resname(
    resids: np.ndarray, resnames: np.ndarray, resids_subset: np.ndarray, resname: str
):
    """Filter the residue IDs by residue name."""
    indices = np.searchsorted(resids, resids_subset)
    result = resids_subset[np.where(resnames[indices] == resname)[0]]
    return set(result)


def filter_lipid_ids_by_resname(
    database, lipid_ids: np.ndarray, lipid_resname: str
) -> np.ndarray:
    """Filter lipid IDs by residue name.

    Parameters
    ----------
    lipid_ids : np.ndarray
        An array of lipid IDs.
    lipid_resname : str
        The residue name to filter by.

    Returns
    -------
    np.ndarray
        An array of filtered lipid IDs.
    """
    sorted_lipid_ids = np.sort(lipid_ids)
    sorted_indices = np.searchsorted(database.residues.resids, sorted_lipid_ids)
    mask = np.zeros(database.residues.resids.shape, dtype=bool)
    mask[sorted_indices] = True
    filtered_resnames = sorted_lipid_ids[
        database.residues.resnames[mask] == lipid_resname
    ]

    return filtered_resnames


def create_lipid_resname_mask(database, lipid_resname):
    """Create a mask for filtering lipid IDs by residue name."""

    return database.residues.resnames == lipid_resname


def filter_resnames_by_lipid_ids_optimized(lipid_resname_mask, lipid_ids, database):
    """Filter lipid IDs by residue name. This is an optimized version of filter_lipid_ids_by_resname, which requires
    the lipid_resname_mask to be precomputed."""
    sorted_lipid_ids = np.sort(lipid_ids)
    sorted_indices = np.searchsorted(database.residues.resids, sorted_lipid_ids)
    mask = np.zeros(database.residues.resids.shape, dtype=bool)
    mask[sorted_indices] = True
    combined_mask = lipid_resname_mask & mask
    filtered_resnames = sorted_lipid_ids[combined_mask[sorted_indices]]

    return filtered_resnames


def contact_frames_to_binary_array(
    contact_frames: Iterable[int], n_frames: int
) -> np.ndarray:
    """Convert a list of contact frames to a binary array.

    Parameters
    ----------
    contact_frames : Iterable[int]
        A list of contact frames.
    n_frames : int
        The number of frames in the trajectory.

    Returns
    -------
    np.ndarray
        A binary array with ones at the indices corresponding to the contact frames.
    """
    binary_array = np.zeros(n_frames)
    binary_array[contact_frames] = 1

    return binary_array


def count_contiguous_segments(arr: np.ndarray) -> np.ndarray:
    """Count the number of contiguous segments of ones in a binary array.

    Parameters
    ----------
    arr : array_like
        A binary array.

    Returns
    -------
    np.ndarray
        An array of segment lengths.
    """
    if np.all(arr == 0):
        return np.array([])

    padded_arr = np.concatenate(([0], arr, [0]))

    start_indices = np.where(np.diff(padded_arr) == 1)[0]
    end_indices = np.where(np.diff(padded_arr) == -1)[0]

    segment_lengths = end_indices - start_indices

    return segment_lengths


def fast_contiguous_segment_lengths(arr, multiplier: float = 1.0) -> np.ndarray:
    """Compute the lengths of contiguous segments of indices in the input array.

    Parameters
    ----------
    arr : Iterable[int]
        A sorted list of indices.

    Returns
    -------
    np.ndarray
        An array of contiguous segment lengths.
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


def index_of_ones(arr: np.ndarray) -> np.ndarray:
    """Return the indices of ones in a binary array.

    Parameters
    ----------
    arr : array_like
        A binary array.

    Returns
    -------
    np.ndarray
        An array of indices.
    """

    return np.where(arr == 1)[0]


def compute_lipid_durations(
    database,
    contact_frames: Dict[int, List[int]],
    lipid_resname: str,
    n_frames: int,
    multiplier: float = 1,
) -> np.ndarray:
    """Compute the duration of lipid contacts. Slower implementation.
    See ContactDurations for a faster implementation.

    Parameters
    ----------
    contact_frames : Iterable[int]
        A list of contact frames.
    n_frames : int
        The number of frames in the trajectory.
    lipid_resname : str
        The residue name of the lipid to compute durations for.

    Returns
    -------
    np.ndarray
        An array of lipid contact durations.
    """

    ids_to_filter = np.array(list(contact_frames.keys()))
    lipid_ids = filter_lipid_ids_by_resname(database, ids_to_filter, lipid_resname)

    durations = [
        contact_frames_to_binary_array(v, n_frames)
        for k, v in contact_frames.items()
        if k in lipid_ids
    ]
    durations = [count_contiguous_segments(v) * multiplier for v in durations]

    return sorted(chain.from_iterable(durations))
