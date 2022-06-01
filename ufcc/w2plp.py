r"""Temporal module to wrap prolintpy
=====================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import numpy as np
from scipy.sparse import hstack
from collections import defaultdict, Counter


class LPContacts(object):
    def __init__(self, mtr_contacts, counts, n_residues_db, frames, PLASMA_LIPIDS, database, timestep, residue='residue'):
        """Create LPContacts class and populate attributes"""
        self.residue = residue
        self.lipids = {}
        self.contacts = {}
        self.occupancy = {}
        for lipid in PLASMA_LIPIDS.keys():
            self.occupancy[lipid] = [0 for frame in range(frames)]
            frames_idx = counts[(counts['# {}'.format(lipid)] != 0) & (counts['ResID'] == self.residue)]['FrameID']
            if len(frames_idx) == 0:                
                self.contacts[lipid] = [0]
                self.lipids[lipid] = []
            else:
                for frame in frames_idx:
                    self.occupancy[lipid][frame] = 1
                global_cont = hstack(mtr_contacts[frames_idx])
                prot_idx, lipid_idx = global_cont.nonzero()
                lipid_idx = lipid_idx[prot_idx == self.residue]
                self.lipids[lipid] = [idx%n_residues_db for idx in lipid_idx]
                self.contacts[lipid] = [x*timestep for x in list((Counter(self.lipids[lipid]).values()))]
                self.lipids[lipid] = np.unique(list(database.selected.residues[self.lipids[lipid]].resindices + 1))
        self.residue += 1 

    def __str__(self):
        return "<ufcc.w2pl2.LPContacts for residue {}>".format(self.residue)


    def __repr__(self):
        return "<ufcc.w2plp.LPContacts for residue {}>".format(self.residue)
