from typing import Iterable, List, Dict
from itertools import chain

import numpy as np

def filter_lipid_ids_by_resname(database, lipid_ids: np.ndarray, lipid_resname: str) -> np.ndarray:
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
    sorted_indices = np.searchsorted(database.selected.residues.resids, sorted_lipid_ids)
    mask = np.zeros(database.selected.residues.resids.shape, dtype=bool)
    mask[sorted_indices] = True
    filtered_resnames = sorted_lipid_ids[database.selected.residues.resnames[mask] == lipid_resname]

    return filtered_resnames


def create_lipid_resname_mask(database, lipid_resname):
    """Create a mask for filtering lipid IDs by residue name. """
    
    return database.selected.residues.resnames == lipid_resname

def filter_resnames_by_lipid_ids_optimized(lipid_resname_mask, lipid_ids, database):
    """Filter lipid IDs by residue name. This is an optimized version of filter_lipid_ids_by_resname, which requires
    the lipid_resname_mask to be precomputed."""
    sorted_lipid_ids = np.sort(lipid_ids)
    sorted_indices = np.searchsorted(database.selected.residues.resids, sorted_lipid_ids)
    mask = np.zeros(database.selected.residues.resids.shape, dtype=bool)
    mask[sorted_indices] = True
    combined_mask = lipid_resname_mask & mask
    filtered_resnames = sorted_lipid_ids[combined_mask[sorted_indices]]

    return filtered_resnames

def contact_frames_to_binary_array(contact_frames: Iterable[int], n_frames: int) -> np.ndarray:
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

def compute_lipid_durations(database, contact_frames: Dict[int, List[int]], lipid_resname: str, n_frames: int, multiplier: float = 1) -> np.ndarray:
    """Compute the duration of lipid contacts.
    
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
    
    durations = [contact_frames_to_binary_array(v, n_frames) for k, v in contact_frames.items() if k in lipid_ids]
    durations = [count_contiguous_segments(v) * multiplier for v in durations]

    return sorted(chain.from_iterable(durations))
