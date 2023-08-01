import numpy as np

from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase

class SerialDistances(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class.
    """
    def __init__(
        self, universe, query, database, lipid_id, residue_id, frame_filter, **kwargs
    ):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.frame_filter = frame_filter
        frame_range = np.arange(len(self.frame_filter))
        self.frame_mapping = {k: v for k, v in zip(self.frame_filter, frame_range)}

        self.lipid_atomgroup = self.database.select_atoms(f"resid {lipid_id}")
        self.resid_atomgroup = self.query.select_atoms(f"resid {residue_id}")
        self.lipid_atomnames = self.lipid_atomgroup.names.tolist()
        self.resid_atomnames = self.resid_atomgroup.names.tolist()
        self.result_array = None
        self.distance_array = None

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

    def _prepare(self):
        self.result_array = np.zeros(
            (
                len(self.frame_filter),
                self.lipid_atomgroup.n_atoms,
                self.resid_atomgroup.n_atoms,
            )
        )

    def _single_frame(self):
        if self._frame_index in self.frame_filter:
            r = distances.distance_array(
                self.lipid_atomgroup.positions,
                self.resid_atomgroup.positions,
                box=self.database.universe.dimensions,
            )
            # print ('frame iterator: ', self._frame_index)
            self.result_array[self.frame_mapping[self._frame_index]] = r

    def _conclude(self):
        self.distance_array = np.mean(self.result_array, axis=0)
        del self.result_array
