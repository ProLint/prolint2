"""Exact contact storage and aggregation.

This module provides the ExactContacts class for storing contacts
with exact frame-level information and computing duration-based metrics.
"""

import logging
from typing import List, Dict, Callable, Union

import numpy as np

from prolint.contacts.base import BaseContactStore
from prolint.utils.utils import fast_contiguous_segment_lengths

logger = logging.getLogger(__name__)


class ExactContacts(BaseContactStore):
    """Exact contact storage with duration-based metric computation.

    Stores contacts at frame-level precision and computes metrics
    based on contiguous contact durations (binding events).

    Parameters
    ----------
    ts : Universe
        MDAnalysis Universe instance.
    contact_frames : dict
        Nested dict mapping residue_id -> database_id -> list of frame indices.
    norm_factor : float, default=1.0
        Normalization factor for duration calculations.

    Attributes
    ----------
    contacts : dict
        Contact durations organized by residue, database type, and molecule ID.
    contact_frames : dict
        Raw frame indices where contacts occur.

    Examples
    --------
    >>> contacts = universe.compute_contacts(cutoff=7.0)
    >>> occupancy = contacts.compute_metric("occupancy", target_resname="CHOL")
    >>> mean_duration = contacts.compute_metric("mean")

    See Also
    --------
    BaseContactStore : Abstract base class
    ComputedContacts : High-level wrapper for contact results
    """

    def run(self, database_resnames: Union[str, List] = None) -> None:
        """Aggregate contact frames into contact durations.

        Processes raw contact frame indices into contiguous binding events
        (durations) for each residue-molecule pair. Results are stored in
        the ``contacts`` attribute.

        Parameters
        ----------
        database_resnames : str or list of str, optional
            Residue names to process. If None, processes all unique
            residue names in the database.
        """
        if database_resnames is None:
            database_resnames = self._universe.database.unique_resnames
        elif isinstance(database_resnames, str):
            database_resnames = [database_resnames]

        logger.debug(
            "Processing %d residues for %d database types",
            len(self.contact_frames),
            len(database_resnames),
        )

        for residue, contact_frame in self.contact_frames.items():
            for database_resname in database_resnames:
                result = self.compute_database_durations(
                    contact_frame, database_resname
                )
                if len(result) > 0:
                    self._contacts[residue][database_resname] = result

        logger.debug(
            "Aggregation complete: %d residues with contacts",
            len(self._contacts),
        )

    def compute_metric(self, metric: str, target_resname=None):
        """Compute a metric across all contacts.

        Parameters
        ----------
        metric : {"max", "sum", "mean", "occupancy"}
            Metric to compute:
            - "occupancy": Fraction of frames with contact
            - "mean": Mean contact duration
            - "max": Maximum contact duration
            - "sum": Total contact duration
        target_resname : str, optional
            Filter by database residue name (e.g., "CHOL").

        Returns
        -------
        dict
            Nested dict with structure:
            {residue_id: {database_name: {"global": value, "per_id": {id: value}}}}
        """
        # Pre-fetch functions and values for speed
        nframes_inv = 1.0 / self._universe.trajectory.n_frames
        is_occupancy = False

        # Get numpy aggregation function once (avoid repeated getattr)
        if metric == "max":
            np_agg_func = np.max
        elif metric == "sum":
            np_agg_func = np.sum
        elif metric == "mean":
            np_agg_func = np.mean
        elif metric == "occupancy":
            is_occupancy = True
        else:
            raise ValueError(f"Unknown metric: {metric}")

        computed_results = {}

        for residue, database_data in self._contacts.items():
            residue_results = {}
            residue_contact_frames = self.contact_frames.get(residue, {})

            for database_name, database_contacts in database_data.items():
                if target_resname is not None and database_name != target_resname:
                    continue

                if not database_contacts:
                    continue

                # Compute per-id values
                if is_occupancy:
                    # Use contact_frames directly (not scaled durations) for correct occupancy
                    per_id_values = {}
                    all_frames = set()
                    for database_id in database_contacts:
                        frames = residue_contact_frames.get(database_id)
                        if frames is not None:
                            per_id_values[database_id] = len(frames) * nframes_inv
                            all_frames.update(frames)
                        else:
                            per_id_values[database_id] = 0.0
                    # Global occupancy: unique frames with at least one contact
                    global_value = len(all_frames) * nframes_inv
                else:
                    # Compute per-id and collect values in one pass
                    per_id_values = {}
                    values_for_global = []
                    for database_id, durations in database_contacts.items():
                        val = float(np_agg_func(durations))
                        per_id_values[database_id] = val
                        values_for_global.append(val)
                    global_value = float(np_agg_func(values_for_global))

                residue_results[database_name] = {
                    "global": global_value,
                    "per_id": per_id_values,
                }

            if residue_results:
                computed_results[residue] = residue_results

        return computed_results

    def apply_function(self, func: Callable, target_resname=None):
        """Apply a custom function to contact duration arrays.

        Parameters
        ----------
        func : callable
            Function that takes an array of durations and returns a value.
        target_resname : str, optional
            Filter by database residue name.

        Returns
        -------
        dict
            Function results organized by residue and database ID.

        Examples
        --------
        >>> # Custom metric: number of binding events
        >>> n_events = contacts.apply_function(len, target_resname="CHOL")
        """
        computed_results = {}
        for residue, database_data in self._contacts.items():
            computed_results[residue] = {}
            for database_name, database_contacts in database_data.items():
                if target_resname is None or database_name == target_resname:
                    computed_contacts_per_id = {
                        database_id: func(contact_array)
                        for database_id, contact_array in database_contacts.items()
                    }
                    computed_results[residue][database_name] = computed_contacts_per_id
        return computed_results

    def compute_database_durations(
        self, contact_frame: Dict[int, List[int]], database_resname: str
    ) -> Dict[int, np.ndarray]:
        """Compute contact durations for a specific database residue type.

        Parameters
        ----------
        contact_frame : dict
            Mapping of database_id -> list of frame indices.
        database_resname : str
            Residue name to filter by.

        Returns
        -------
        dict
            Mapping of database_id -> array of contact durations.
        """
        ids_to_filter = np.array(list(contact_frame.keys()))
        database_ids = set(
            self._universe.database.filter_resids_by_resname(
                ids_to_filter, database_resname
            )
        )

        durations = {}
        for k, arr in contact_frame.items():
            if k in database_ids:
                durations[k] = fast_contiguous_segment_lengths(arr, self.norm_factor)

        return durations
