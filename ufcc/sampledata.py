r"""UFCC sample data
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import pathlib

class GIRKDataSample:
    def __init__(self):
        self.path = pathlib.Path('data').absolute()
        self.coordinates = self.path.joinpath('coordinates.gro').absolute()
        self.trajectory = self.path.joinpath('trajectory.xtc').absolute()

GIRK = GIRKDataSample()
