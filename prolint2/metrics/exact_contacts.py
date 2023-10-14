r""":mod:`prolint2.metrics.exact_contacts`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from typing import List, Dict, Callable, Union

from collections import defaultdict

import numpy as np

from prolint2.metrics.base import BaseContactStore
from prolint2.metrics.utils import (
    fast_filter_resids_by_resname,
    fast_contiguous_segment_lengths,
)


class ExactContacts(BaseContactStore):
    """Compute the duration of lipid contacts. This class is used to compute the duration of lipid contacts."""

    def run(self, lipid_resnames: Union[str, List] = None) -> Dict[str, np.ndarray]:
        """Compute the duration of lipid contacts for all lipid types.

        Parameters
        ----------
        lipid_resnames : str, optional
            A list of lipid residue names to compute durations for. If None, durations will be computed for all lipid types.

        Returns
        -------
        Dict[str, np.ndarray]
            A dictionary of lipid contact durations for all lipid types.
            The output is stored in the `self._contacts` attribute.
        """
        if lipid_resnames is None:
            lipid_resnames = self._database_unique_resnames
        elif isinstance(lipid_resnames, str):
            lipid_resnames = [lipid_resnames]

        for residue, contact_frame in self.contact_frames.items():
            for lipid_resname in lipid_resnames:
                result = self.compute_lipid_durations(contact_frame, lipid_resname)
                if len(result) > 0:
                    self._contacts[residue][lipid_resname] = result

    def pooled_results(self, target_lipid_name=None):
        """Pool results for all lipids.

        Parameters
        ----------
        target_lipid_name : str, optional
            The name of the lipid to compute pooled results for. If None, pooled results will be computed for all lipids.

        Returns
        -------
        Dict[str, Dict[str, List[float]]]
            A dictionary of pooled results for all lipids.
        """
        pooled_results = defaultdict(lambda: defaultdict(list))
        for residue, lipid_data in self._contacts.items():
            for lipid_name, lipid_contacts in lipid_data.items():
                if target_lipid_name is None or lipid_name == target_lipid_name:
                    pooled_contact_array = []
                    for lipid_id_contacts in lipid_contacts.values():
                        pooled_contact_array.extend(lipid_id_contacts)
                    pooled_results[residue][lipid_name].extend(pooled_contact_array)
        return pooled_results

    def compute_metric(self, metric: str, target_lipid_name=None):
        """Compute a pre-defined metric for all lipids or a specific lipid.

        Parameters
        ----------
        metric : str
            The metric to compute. Must be one of 'max', 'sum', 'mean'.
        target_lipid_name : str, optional
            The name of the lipid to compute the metric for. If None, the metric will be computed for all lipids.

        Returns
        -------
        Dict[str, Dict[str, Dict[int, float]]]
            A dictionary of computed metrics for all lipids.
        """

        computed_results = defaultdict(lambda: defaultdict(dict))
        for residue, lipid_data in self._contacts.items():
            # computed_results[residue] = {}
            for lipid_name, lipid_contacts in lipid_data.items():
                if target_lipid_name is None or lipid_name == target_lipid_name:
                    computed_contacts_per_id = {
                        lipid_id: getattr(np, metric)(contact_array)
                        for lipid_id, contact_array in lipid_contacts.items()
                    }
                    computed_results[residue][lipid_name] = computed_contacts_per_id
        return computed_results

    def apply_function(self, func: Callable, target_lipid_name=None):
        """Apply a function to all lipids or a specific lipid.

        Parameters
        ----------
        func : Callable
            The function to apply to the lipid contact durations.
        target_lipid_name : str, optional
            The name of the lipid to apply the function to. If None, the function will be applied to all lipids.

        Returns
        -------
        Dict[str, Dict[str, Dict[int, float]]]
            A dictionary of computed metrics for all lipids.

        Example
        -------
        >>> cd = ExactContacts(...)
        >>> cd.run()
        >>> cd.apply_function(np.mean)
        >>> cd.apply_function(np.max, target_lipid_name='DOPC')
        >>> cd.apply_function(lambda x: np.mean(x) / np.max(x), target_lipid_name='DOPC')
        """
        computed_results = {}
        for residue, lipid_data in self._contacts.items():
            computed_results[residue] = {}
            for lipid_name, lipid_contacts in lipid_data.items():
                if target_lipid_name is None or lipid_name == target_lipid_name:
                    computed_contacts_per_id = {
                        lipid_id: func(contact_array)
                        for lipid_id, contact_array in lipid_contacts.items()
                    }
                    computed_results[residue][lipid_name] = computed_contacts_per_id
        return computed_results

    def compute_lipid_durations(
        self, contact_frame: Dict[int, List[int]], lipid_resname: str
    ) -> np.ndarray:
        """Compute the duration of lipid contacts.

        Parameters
        ----------
        contact_frame : Dict[int, List[int]]
            A dictionary of contact frames.
        lipid_resname : str
            The residue name of the lipid to compute durations for.

        Returns
        -------
        np.ndarray
            An array of lipid contact durations.
        """

        ids_to_filter = np.array(list(contact_frame.keys()))
        lipid_ids = fast_filter_resids_by_resname(
            self._resids, self._resnames, ids_to_filter, lipid_resname
        )

        durations = {}
        for k, arr in contact_frame.items():
            if k in lipid_ids:
                durations[k] = fast_contiguous_segment_lengths(arr, self.norm_factor)

        return durations
