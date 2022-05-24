"""
Ultr-Fast Contacts Calculation (UFCC)
=======================================

:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License

UFCC calculates de distance-based contacts between two references.
"""

import numpy as np
import warnings

from MDAnalysis.analysis.base import AnalysisBase
from .nsgrid import FastNS

# import logging
# MDAnalysis.start_logging()

# logger = logging.getLogger("MDAnalysis.MDAKit.membrane_curvature")


class Contacts(AnalysisBase):
    r"""
    Class to get the distance-based contacts starting from two AtomGroups.
    """

    def __init__(self, universe, query, database, 
                radius, wrap=True, **kwargs):

        super().__init__(universe.universe.trajectory, **kwargs)
        self.query = query
        self.database = database
        self.radius = radius

        # Raise if selection doesn't exist
        if len(self.query) == 0 or len(self.database) == 0:
            raise ValueError("Invalid selection. Empty AtomGroup(s).")

        # # Apply wrapping coordinates
        # if not self.wrap:
        #     # Warning
        #     msg = (" `wrap == False` may result in inaccurate calculation "
        #            "of membrane curvature. Surfaces will be derived from "
        #            "a reduced number of atoms. \n "
        #            " Ignore this warning if your trajectory has "
        #            " rotational/translational fit rotations! ")
        #     warnings.warn(msg)
        #     logger.warn(msg)

    def _prepare(self):
        # Initialize empty np.array with results
        self.results.contacts = {}

    def _single_frame(self):
        # Get the results and populate the results dictionary
        gridsearch = FastNS(self.radius, self.database.positions, box=self.database.dimensions, pbc=True)
        result = gridsearch.search(self.query.positions)
        pairs = result.get_neighbours()
        contacts_per_idx = []
        flag = 0
        for pair_idx in range(len(pairs)):
            if pair_idx == 0:
                contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
            else:
                if pairs[pair_idx][0] == pairs[pair_idx-1][0]:
                    contacts_per_idx[flag].append(pairs[pair_idx][1])
                else:
                    flag += 1
                    contacts_per_idx.append([pairs[pair_idx][0], pairs[pair_idx][1]])
        self.results.contacts[self._frame_index] = contacts_per_idx


