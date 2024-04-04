r""":mod:`prolint2.computers.contacts`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

from collections import defaultdict
import numpy as np

from MDAnalysis.lib.nsgrid import FastNS

from prolint2.computers.base import ContactComputerBase
from prolint2.utils.utils import fast_unique_comparison

import os
import configparser

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini"))
parameters_config = config["Parameters"]


class SerialContacts(ContactComputerBase):
    """
    Class to get the distance-based contacts starting from two AtomGroups
    using a *serial* approach.

    It inherits from the MDAnalysis AnalysisBase class.

    :param universe: The MDAnalysis Universe object.
    :type universe: MDAnalysis.Universe
    :param query: AtomGroup representing the query.
    :type query: MDAnalysis.AtomGroup
    :param database: AtomGroup representing the database.
    :type database: MDAnalysis.AtomGroup
    :param cutoff: The distance cutoff (default from parameters_config).
    :type cutoff: float
    :param kwargs: Additional keyword arguments.
    :type kwargs: dict
    """

    def __init__(
        self,
        universe,
        query,
        database,
        cutoff=float(parameters_config["cutoff"]),
        **kwargs
    ):
        """
        Initialize the SerialContacts object.

        :param universe: The MDAnalysis Universe object.
        :param query: AtomGroup representing the query.
        :param database: AtomGroup representing the database.
        :param cutoff: The distance cutoff (default from parameters_config).
        :param kwargs: Additional keyword arguments.

        This method initializes the SerialContacts object. It sets up the query, database, and cutoff, and performs input validation.

        :raises ValueError: If the input selection is empty or the cutoff is not greater than 0.
        """
        super().__init__(universe.universe.trajectory, **kwargs)

        self.query = query
        self.database = database
        self.cutoff = cutoff

        self.q_resids = self.query.resids
        self.db_resids = self.database.resids
        self.db_resnames = self.database.resnames

        self.contacts = None
        self.contact_frames = defaultdict(lambda: defaultdict(list))

        self._validate_inputs()

    def _validate_inputs(self):
        """
        Validate the inputs.

        This method validates the inputs by checking if the query and database selections are not empty and if the cutoff is greater than 0.

        :raises ValueError: If the input selection is empty or the cutoff is not greater than 0.
        """
        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        if self.cutoff <= 0:
            raise ValueError("The cutoff must be greater than 0.")

    def _get_residue_lipid_info(self, pair):
        """
        Get the residue and lipid information for a given pair.

        :param pair: A pair of indices.

        :return: A tuple containing the residue ID, lipid ID, and lipid name.
        :rtype: tuple
        """
        residue_id = self.q_resids[pair[0]]
        lipid_id = self.db_resids[pair[1]]
        lipid_name = self.db_resnames[pair[1]]
        return residue_id, lipid_id, lipid_name

    def _compute_pairs(self):
        """
        Compute the pairs of residues and lipids that are within the cutoff distance.

        :return: An array of pairs.
        :rtype: numpy.ndarray
        """

        if self.database.dimensions is None:
            dim_x = np.max(self.database.universe.atoms.positions[:, 0]) - np.min(self.database.universe.atoms.positions[:, 0])
            dim_y = np.max(self.database.universe.atoms.positions[:, 1]) - np.min(self.database.universe.atoms.positions[:, 1])
            dim_z = np.max(self.database.universe.atoms.positions[:, 2]) - np.min(self.database.universe.atoms.positions[:, 2])
            self.database.dimensions = np.array([dim_x, dim_y, dim_z, 90, 90, 90])
        
        gridsearch = FastNS(
            self.cutoff, self.database.positions, box=self.database.dimensions, pbc=True
        )
        result = gridsearch.search(self.query.positions)
        pairs = result.get_pairs()

        return pairs

    def _single_frame(self):
        """
        Compute the contacts for a single frame.

        This method computes contacts for a single frame by iterating through pairs of residues and lipids within the cutoff distance.

        It updates the `contact_frames` attribute.
        """
        pairs = self._compute_pairs()

        q_resid_indices = pairs[:, 0]
        db_resid_indices = pairs[:, 1]
        residue_ids = self.q_resids[q_resid_indices]
        lipid_ids = self.db_resids[db_resid_indices]
        lipid_names = self.db_resnames[db_resid_indices]

        residue_ids, lipid_ids, lipid_names = fast_unique_comparison(
            residue_ids, lipid_ids, lipid_names
        )

        existing_pairs = set()
        for unique_data in zip(residue_ids, lipid_ids, lipid_names):
            residue_id, lipid_id, _ = unique_data
            if (residue_id, lipid_id) not in existing_pairs:
                existing_pairs.add((residue_id, lipid_id))
                self.contact_frames[residue_id][lipid_id].append(self._frame_index)

    # def _conclude(self):
    #     contacts = ExactContacts(self.query.universe, self.contact_frames)
    #     contacts.run()
    #     self.contacts = contacts
