r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os

# to get the paths relative to the root of the package
from prolint2 import get_data


class GIRKDataSample:
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the GIRK channel using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import GIRK

        target_system = PL2(GIRK.coordinates, GIRK.trajectory)

    """

    def __init__(self):
        self.path = os.path.abspath(get_data())
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")


# initializing GIRK sample data
GIRK = GIRKDataSample()
