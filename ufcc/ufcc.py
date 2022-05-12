# Daniel P. Ramirez & Besian I. Sejdiu
# Prolint: A tool to analyze and visualize lipid-protein interactions.
#

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import _ResidueStringAttr



def get_atom_group(selection):
    """
    Cast an MDAnalysis.Atom, MDAnalysis.Residue, or MDAnalysis.ResidueGroup to AtomGroup.
    Parameters
    ----------
    selection: MDAnalysis.Atom, MDAnalysis.Residue or MDAnalysis.ResidueGroup
        atoms to cast
    Returns
    -------
    MDAnalysis.AtomGroup
    """
    assert isinstance(
        selection,
        (
            mda.core.groups.Residue,
            mda.core.groups.ResidueGroup,
            mda.core.groups.Atom,
            mda.core.groups.AtomGroup,
        ),
    ), "central_species must be one of the preceding types"
    u = selection.universe
    if isinstance(selection, (mda.core.groups.Residue, mda.core.groups.ResidueGroup)):
        selection = selection.atoms
    if isinstance(selection, mda.core.groups.Atom):
        selection = u.select_atoms(f"index {selection.index}")
    return selection


class MacrosClass(_ResidueStringAttr):
    attrname = 'macros'
    singular = 'macro'

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments):
        return np.array(['other'] * n_residues, dtype=object)


class UFCC(object):
    """Base class for getting topology information. It reads an MDAnalysis Universe
    and extracts useful information. 

    Attributes
    ----------
    universe : MDAnalysis universe with all different kind of topologies
    """
    def __init__(self, structure, trajectory):
        self.universe = mda.Universe(structure, trajectory)
        self.atoms = self.universe.atoms
        self.list_atoms = np.unique(self.universe.atoms)
        self.residues = self.universe.residues
        self.list_residues = np.unique(self.universe.residues)
        # self.molecules = self.universe.atoms.moltypes ## this information is only on the tpr, not on the gro
        self.universe.add_TopologyAttr('macros')
        self.universe.select_atoms('protein').residues.macros = 'protein'
        self.universe.select_atoms('not protein').residues.macros = 'membrane'
        self.macros = np.unique(self.universe.residues.macros)

    def get_atom_group(selection):
        """
        Cast an MDAnalysis.Atom, MDAnalysis.Residue, or MDAnalysis.ResidueGroup to AtomGroup.
        Parameters
        ----------
        selection: MDAnalysis.Atom, MDAnalysis.Residue or MDAnalysis.ResidueGroup
            atoms to cast
        Returns
        -------
        MDAnalysis.AtomGroup
        """
        assert isinstance(
            selection,
            (
                mda.core.groups.Residue,
                mda.core.groups.ResidueGroup,
                mda.core.groups.Atom,
                mda.core.groups.AtomGroup,
            ),
        ), "central_species must be one of the preceding types"
        u = selection.universe
        if isinstance(selection, (mda.core.groups.Residue, mda.core.groups.ResidueGroup)):
            selection = selection.atoms
        if isinstance(selection, mda.core.groups.Atom):
            selection = u.select_atoms(f"index {selection.index}")
        return selection  



