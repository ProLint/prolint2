from typing import List, Dict, Literal

from collections import defaultdict
from itertools import chain

import numpy as np

from prolint2.metrics.utils import fast_filter_resids_by_resname, fast_contiguous_segment_lengths


class ContactDurations:
    """Compute the duration of lipid contacts. This class is used to compute the duration of lipid contacts. """
    def __init__(self, ts, contacts, custom_multiplier: float = None, unit: Literal['us', 'ns'] = 'us'):
        self.n_frames = ts.trajectory.n_frames
        self.unit = unit
        self.unit_divisor = 1000000 if unit == 'us' else 1000
        self.multiplier = custom_multiplier if custom_multiplier is not None else ts.trajectory.dt / self.unit_divisor
        self.contact_frames = contacts.contact_frames

        # stored for performance reasons. MDAnalysis `resids` and `resnames` accesses are slow.
        self._resids = ts.database.residues.resids
        self._resnames = ts.database.residues.resnames
        self._lipid_resname_mask = {k: self._create_lipid_resname_mask(ts.database, k) for k in ts.database.unique_resnames}
        self._durations = None

    def compute(self, contact_frame: Dict[int, List[int]], lipid_resname: str) -> np.ndarray:
        """Compute the duration of lipid contacts.

        Parameters
        ----------
        contact_frames : Dict[int, List[int]]
            A dictionary of contact frames.
        lipid_resname : str
            The residue name of the lipid to compute durations for.

        Returns
        -------
        np.ndarray
            An array of lipid contact durations.
        """

        return self.compute_lipid_durations(contact_frame, lipid_resname)

    def compute_all(self, lipid_resnames: str = None) -> Dict[str, np.ndarray]:
        """Compute the duration of lipid contacts for all lipid types.

        Parameters
        ----------
        lipid_resnames : str, optional
            A list of lipid residue names to compute durations for. If None, durations will be computed for all lipid types.

        Returns
        -------
        Dict[str, np.ndarray]
            A dictionary of lipid contact durations for all lipid types.
        """
        if lipid_resnames is None:
            lipid_resnames = self._lipid_resname_mask.keys()
        elif isinstance(lipid_resnames, str):
            lipid_resnames = [lipid_resnames]

        results = defaultdict(dict)
        for residue, contact_frame in self.contact_frames.items():
            for lipid_resname in lipid_resnames:
                # results[lipid_resname][residue] = self.compute(contact_frame, lipid_resname)
                result = self.compute(contact_frame, lipid_resname)
                if len(result) > 0:
                    results[lipid_resname][residue] = result

        self._durations = results

    def compute_lipid_durations(self, contact_frames: Dict[int, List[int]], lipid_resname: str) -> np.ndarray:
        """Compute the duration of lipid contacts.
        
        Parameters
        ----------
        contact_frames : Dict[int, List[int]]
            A dictionary of contact frames.
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
        lipid_ids = fast_filter_resids_by_resname(self._resids, self._resnames, ids_to_filter, lipid_resname)
        
        durations = []
        for k, arr in contact_frames.items():
            if k in lipid_ids:
                durations.append(fast_contiguous_segment_lengths(arr, self.multiplier))
        
        return sorted(chain.from_iterable(durations))

    @staticmethod
    def _create_lipid_resname_mask(database, lipid_resname):
        # set the lipid_resname_mask to True for all lipids with the given resname
        return database.residues.resnames == lipid_resname

    @property
    def results(self):
        """Get the computed durations. """
        if self._durations is None:
            raise ValueError('No durations have been computed yet. Call compute_all() first.')
        return self._durations
