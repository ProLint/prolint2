# Daniel P. Ramirez & Besian I. Sejdiu
# Prolint: A tool to analyze and visualize lipid-protein interactions.
#

import MDAnalysis as mda

class Contacts(object):
    """Base class for getting the contacts. 

    Atributes
    ----------
    pairs : (idx_query, idx_haystack) pairs.
    """
    def __init__(self):
        self.pairs = 'TODO'