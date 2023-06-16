r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import os

# to get the paths relative to the root of the package

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
        self.path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/GIRK/")
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


# initializing GIRK sample data
GIRK = GIRKDataSample()


class COX1DataSample:
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the COX1 enzyme using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import COX1

        target_system = PL2(COX1.coordinates, COX1.trajectory)

    """

    def __init__(self):
        self.path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/COX1/")
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


# initializing COX1 sample data
COX1 = COX1DataSample()


class SMODataSample:
    """
    Class to add sample data in order to test the installation of the package

    It uses the coordinates and trajectory of a simulation of the SMO receptor using the MARTINI force field.

    Example
    -------
    Import the prolint2 library and use the sample data as follows::

        from prolint2 import PL2
        from prolint2.sampledata import SMO

        target_system = PL2(SMO.coordinates, SMO.trajectory)

    """

    def __init__(self):
        self.path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data/SMO/")
        self.coordinates = os.path.join(self.path, "coordinates.gro")
        self.trajectory = os.path.join(self.path, "trajectory.xtc")
        self.contacts = os.path.join(self.path, "contacts.csv")


# initializing SMO sample data
SMO = SMODataSample()