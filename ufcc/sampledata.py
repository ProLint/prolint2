r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os

from ufcc import get_data


class GIRKDataSample:
    def __init__(self):
        self.path = os.path.abspath(get_data())
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")


GIRK = GIRKDataSample()
