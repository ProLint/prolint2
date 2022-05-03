# Daniel P. Ramirez & Besian I. Sejdiu
# Prolint: A tool to analyze and visualize lipid-protein interactions.
#

import MDAnalysis as mda

class UFCC(object):
    """Base class for getting topology information. It reads an MDAnalysis Universe
    and extracts useful information. 

    Attributes
    ----------
    top : MDAnalysis universe with all different kind of topologies
    """
    def __init__(self, structure, trajectory):
        self.top = mda.Universe(structure, trajectory)


