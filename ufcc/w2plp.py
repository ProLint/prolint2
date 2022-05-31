r"""Temporal module to wrap prolintpy
=====================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import numpy as np

class LPContacts(object):
    def __init__(self, mtr_contacts, counts, PLASMA_LIPIDS, database, timestep, residue='residue'):
        """Create LPContacts class and populate attributes"""
        self.residue = residue
        self.lipids = {}
        self.contacts = {}
        self.occupancy = {}
        frames = database.selected.universe.trajectory.n_frames
        for lipid in PLASMA_LIPIDS.keys():
            self.contacts[lipid] = []
            self.lipids[lipid] = []
            self.occupancy[lipid] = [0 for frame in range(frames)]
            mask = counts['# {}'.format(lipid)] != 0 
            temp = counts[mask]
            mask = temp['ResID'] == self.residue
            frames_idx = temp[mask]['FrameID']
            if len(frames_idx) == 0:
                pass
                self.contacts[lipid] = self.lipids[lipid] + [0]
            else:
                for frame in frames_idx:
                    self.occupancy[lipid][frame] = 1
                    lipid_idx = mtr_contacts[frame].nonzero()[1][mtr_contacts[frame].nonzero()[0] == self.residue]
                    self.lipids[lipid] = self.lipids[lipid] + [database.selected.residues[idx].resindex for idx in lipid_idx]
                    self.contacts[lipid] = self.contacts[lipid] + [timestep for x in range(len(lipid_idx))]

    def __str__(self):
        return "<prolintpy.LPContacts for residue {}>".format(self.residue)


    def __repr__(self):
        return "<prolintpy.LPContacts for residue {}>".format(self.residue)