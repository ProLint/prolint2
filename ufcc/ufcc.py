r"""UFCC base class --- :mod:`ufcc.UFCC`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License

UFCC calculates de distance-based contacts between two references.

The class and its methods
-------------------------
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import ResidueStringAttr
from .contacts import Contacts


class MacrosClass(ResidueStringAttr):
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

    def __init__(self, structure, trajectory, add_lipid_types = []):
        self.atoms = mda.Universe(structure, trajectory).atoms
        self.residues = self.atoms.residues
        self.atoms.universe.add_TopologyAttr('macros')

        lipid_types = ['POPC', 'DPPC', 'DOPC', 'CHOL', 'CHL1', 'POPS', 'POPE']
        lipid_types = lipid_types + add_lipid_types
        not_protein_restypes = np.unique(self.atoms.select_atoms('not protein').residues.resnames)
        membrane_restypes = []
        for type in lipid_types:
            if type in not_protein_restypes:
                membrane_restypes.append('resname ' + type)
        if len(membrane_restypes) == 1:
            membrane_sel = membrane_restypes[0]
        elif len(membrane_restypes) > 1:
            membrane_sel = membrane_restypes[0]
            for type in membrane_restypes[1:]:
                membrane_sel = membrane_sel + ' or ' + type
        else:
            print('There are not lipid residues in your system')

        protein_sel = 'protein'
        if len(self.atoms.select_atoms(protein_sel).segments) > 1 and self.atoms.select_atoms(protein_sel).segments.n_atoms == self.atoms.select_atoms(protein_sel).n_atoms:
            for segment_idx in range(len(self.atoms.select_atoms(protein_sel).segments)):
                self.atoms.select_atoms(protein_sel).segments[segment_idx].residues.macros = 'protein' + str(segment_idx)
        else:
            # Get start and end indices of proteins in the system.
            # The assumption here is that proteins are ordered and the start residue of the next
            # protein is always smaller than the last residue of the previous protein.
            resseq = self.atoms.select_atoms(protein_sel).residues.resindices
            p0 = resseq[0]
            # system_proteins contains the start and end indices of all proteins.
            fi_li = []
            fi = 0
            for li, p in enumerate(resseq):
                if p < p0:
                    fi_li.append((fi, li-1))
                    fi = li
                p0 = p
            fi_li.append((fi, li))
            
            for idx, values in enumerate(fi_li):
                # first and last index
                fi = values[0]
                li = values[1]
                self.atoms.select_atoms(protein_sel).residues[list(range(fi, li+1))].residues.macros = 'protein' + str(idx)

        # TODO
        # Add merge chains and options to change the name of the proteins.

        self.atoms.select_atoms(membrane_sel).residues.macros = 'membrane'
        self.list_macros = list(np.unique(self.atoms.residues.macros))
        self.query = QueryProteins(self.atoms.select_atoms(protein_sel))
        self.database = MembraneDatabase(self.atoms.select_atoms(membrane_sel))
        self.contacts = Contacts(self.query, self.database)

    def __str__(self):
        return "Base class for the contacts routines."

    def __repr__(self):
        return "Base class for the contacts routines."


class BasicGroup(object):

    def __init__(self, whole):
        self.AG = whole
        self.whole = whole


    def select(self, selection='all'):
        """
        Cast an MDAnalysis.Atom, MDAnalysis.Residue, or MDAnalysis.ResidueGroup, or str syntax to AtomGroup.
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
            (str, np.ndarray, mda.core.groups.Residue, mda.core.groups.ResidueGroup, mda.core.groups.Atom,
             mda.core.groups.AtomGroup),
        ), "the selection must be one of the preceding types"
        if isinstance(selection, (mda.core.groups.Residue, mda.core.groups.ResidueGroup)):
            selection = selection.atoms
        elif isinstance(selection, mda.core.groups.Atom):
            selection = self.whole.select_atoms(f"index {selection.index}")
        elif isinstance(selection, np.ndarray):
            selection = self.whole.atoms[selection]
        elif isinstance(selection, str):
            selection = self.whole.atoms.select_atoms(selection)
        self.AG = selection


class MembraneDatabase(BasicGroup):
    """Class for the membrane group."""

    def __init__(self, whole):
        super().__init__(whole)

    def lipid_types(self):
        """Get the names of all lipids that will be analyzed.
        Returns
        -------
        names : array of lipid names
        """
        if not isinstance(self.AG, mda.core.groups.AtomGroup):
            return np.array([])
        else:
            return np.unique(self.AG.resnames)

    def lipid_count(self):
        """Get the name and count of each lipid that will be analyzed.
        Returns
        -------
        lipic_count : dictionary
            key:value corresponds to lipid_name:count.
        """
        lc = {}
        lipids = self.lipid_types()
        for lipid in lipids:
            lc[lipid] = len(self.AG.residues[self.AG.residues.resnames == lipid])
        return lc

    def __str__(self):
        if not isinstance(self.AG, mda.core.groups.AtomGroup):
            return "<prolintpy.MembraneDatabase containing 0 atoms>"
        else:
            return "<prolintpy.MembraneDatabase containing {} atoms>".format(self.AG.n_atoms)

    def __repr__(self):
        if not isinstance(self.AG, mda.core.groups.AtomGroup):
            return "<prolintpy.MembraneDatabase containing 0 atoms>"
        else:
            return "<prolintpy.MembraneDatabase containing {} atoms>".format(self.AG.n_atoms)


class QueryProteins(BasicGroup):
    """
    Proteins stores information about the protein composition of the system.
    Gromacs coordinate files do not contain information on protein names and their
    count. This class tries to calculate and extract this information from the topology.
    The class has a 'system_proteins' method that can be used to retrieve information
    about each protein specifically.
    """

    def __init__(self, whole):
        """Create a new ProLint Proteins class."""
        super().__init__(whole)

    def list_proteins(self):
        return np.unique(self.whole.residues.macros)

    def __str__(self):
        if not isinstance(self.AG, mda.core.groups.AtomGroup):
            return "<prolintpy.QueryProteins containing 0 atoms>"
        else:
            return "<prolintpy.QueryProteins containing {} atoms>".format(self.AG.n_atoms)

    def __repr__(self):
        if not isinstance(self.AG, mda.core.groups.AtomGroup):
            return "<prolintpy.QueryProteins containing 0 atoms>"
        else:
            return "<prolintpy.QueryProteins containing {} atoms>".format(self.AG.n_atoms)
