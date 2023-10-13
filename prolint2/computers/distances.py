r""":mod:`prolint2.computers.distances`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import numpy as np

from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase


class SerialDistances(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class.

    :param universe: The MDAnalysis Universe.
    :type universe: MDAnalysis.Universe
    :param query: AtomGroup to query.
    :type query: MDAnalysis.AtomGroup
    :param database: AtomGroup to use as the database.
    :type database: MDAnalysis.AtomGroup
    :param lipid_id: Residue ID for the lipid group.
    :type lipid_id: int
    :param residue_id: Residue ID for the query group.
    :type residue_id: int
    :param frame_filter: List of frame indices to analyze.
    :type frame_filter: list
    :param kwargs: Additional keyword arguments.
    :type kwargs: dict
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
        """
        Prepare data structures for storing distance calculations.

        This method initializes the result_array to store distance results.

        :return: None
        """
        self.result_array = np.zeros(
            (
                len(self.frame_filter),
                self.lipid_atomgroup.n_atoms,
                self.resid_atomgroup.n_atoms,
            )
        )

    def _single_frame(self):
        """
        Perform distance calculations for a single frame.

        This method calculates distances between atoms of lipid and residue groups
        for the current frame and stores the results in result_array.

        :return: None
        """
        if self._frame_index in self.frame_filter:
            r = distances.distance_array(
                self.lipid_atomgroup.positions,
                self.resid_atomgroup.positions,
                box=self.database.universe.dimensions,
            )
            # print ('frame iterator: ', self._frame_index)
            self.result_array[self.frame_mapping[self._frame_index]] = r

    def _conclude(self):
        """
        Conclude the distance analysis and compute the mean distance array.

        This method calculates the mean distance array based on the results stored
        in result_array and deletes the result_array to free up memory.

        :return: None
        """
        self.distance_array = np.mean(self.result_array, axis=0)
        del self.result_array


class TwoPointDistances(AnalysisBase):
    """
    Initialize the TwoPointDistances analysis.

    :param universe: The molecular dynamics universe.
    :type universe: MDAnalysis.Universe
    :param query: The query object for comparison.
    :type query: MDAnalysis.AtomGroup
    :param database: The database object for comparison.
    :type database: MDAnalysis.AtomGroup
    :param lipid_id: The ID of the lipid to consider.
    :type lipid_id: int
    :param residue_id: The ID of the residue to consider.
    :type residue_id: int
    :param lipid_sel: The optional selection for the lipid atom.
    :type lipid_sel: str
    :param residue_sel: The optional selection for the residue atom.
    :type residue_sel: str
    :param unit: The unit for calculating distances ("frame" or "time").
    :type unit: str
    :param **kwargs: Additional keyword arguments.
    :type kwargs: dict

    This method sets up the analysis and initializes necessary attributes.
    """

    def __init__(
        self,
        universe,
        query,
        database,
        lipid_id,
        residue_id,
        lipid_sel=None,
        residue_sel=None,
        unit="frame",
        **kwargs,
    ):
        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.unit = unit
        if lipid_sel is not None and residue_sel is not None:
            self.pointA = self.database.select_atoms(
                f"resid {lipid_id} and name {lipid_sel}"
            )
            self.pointB = self.query.select_atoms(
                f"resid {residue_id} and name {residue_sel}"
            )
        elif lipid_sel is None and residue_sel is None:
            self.pointA = self.database.select_atoms(
                f"resid {lipid_id}"
            ).center_of_mass(compound="residues")
            self.pointB = self.query.select_atoms(f"resid {residue_id}").center_of_mass(
                compound="residues"
            )
        elif lipid_sel is None:
            self.pointA = self.database.select_atoms(
                f"resid {lipid_id}"
            ).center_of_mass(compound="residues")
            self.pointB = self.query.select_atoms(
                f"resid {residue_id} and name {residue_sel}"
            )
        elif residue_sel is None:
            self.pointA = self.database.select_atoms(
                f"resid {lipid_id} and name {lipid_sel}"
            )
            self.pointB = self.query.select_atoms(f"resid {residue_id}").center_of_mass(
                compound="residues"
            )
        self.result_array = None
        self.time_array = None

    def _prepare(self):
        """
        Prepare the analysis.

        This method is called before the analysis starts and prepares the result array.
        """
        self.result_array = np.zeros(self.database.universe.trajectory.n_frames)

    def _single_frame(self):
        """
        Calculate distances for a single frame.

        This method is called for each frame in the trajectory to calculate distances.
        """
        if isinstance(self.pointA, np.ndarray) and isinstance(self.pointB, np.ndarray):
            dist = distances.distance_array(
                self.pointA,
                self.pointB,
                box=self.database.universe.dimensions,
            )
        elif isinstance(self.pointA, np.ndarray):
            dist = distances.distance_array(
                self.pointA,
                self.pointB.positions,
                box=self.database.universe.dimensions,
            )
        elif isinstance(self.pointB, np.ndarray):
            dist = distances.distance_array(
                self.pointA.positions,
                self.pointB,
                box=self.database.universe.dimensions,
            )
        else:
            dist = distances.distance_array(
                self.pointA.positions,
                self.pointB.positions,
                box=self.database.universe.dimensions,
            )
        self.result_array[self._frame_index] = float(dist)

    def _conclude(self):
        """
        Conclude the analysis.

        This method is called after the analysis is complete and sets the time_array attribute.
        """
        if self.unit == "frame":
            self.time_array = self.frames
        elif self.unit == "time":
            self.time_array = self.times
