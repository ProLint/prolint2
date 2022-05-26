r"""UFCC base class --- :mod:`ufcc.UFCC`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License

UFCC calculates de distance-based contacts between two references.

The class and its methods
-------------------------
.. autoclass:: UFCC
    :members:
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

    def __init__(self, structure, trajectory):
        self.atoms = mda.Universe(structure, trajectory).atoms
        self.residues = self.atoms.residues
        self.atoms.universe.add_TopologyAttr('macros')
        self.atoms.select_atoms('protein').residues.macros = 'protein'
        self.atoms.select_atoms('not protein').residues.macros = 'membrane'
        self.list_macros = list(np.unique(self.atoms.residues.macros))
        self.query = QueryProteins(self.atoms.universe)
        self.database = MembraneDatabase(self.atoms.universe)
        self.contacts = Contacts(self.query, self.database)
   
    def __str__(self):
        return "Base class for the contacts routines."

    def __repr__(self):
        return "Base class for the contacts routines."


class BasicGroup(object):
    def __init__(self, universe):
        self.AG = None
        self.universe = universe

    def select(self, selection='all', add_filter = 'all'):
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
        assert isinstance(
            add_filter,
            (str),
        ), "the filter must be always a string"
        if isinstance(selection, (mda.core.groups.Residue, mda.core.groups.ResidueGroup)):
            selection = selection.atoms
        elif isinstance(selection, mda.core.groups.Atom):
            selection = self.universe.select_atoms(f"index {selection.index}")
        elif isinstance(selection, np.ndarray):
            selection = self.universe.atoms[selection]
        elif isinstance(selection, str):
            selection = self.universe.atoms.select_atoms(selection)
        self.AG = selection.select_atoms(add_filter)


class MembraneDatabase(BasicGroup):
    """Class for the membrane group."""
    def __init__(self, universe):
        super().__init__(universe)
    
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

    def __init__(self, universe):
        """Create a new ProLint Proteins class."""
        super().__init__(universe)

    # TODO TODO TODO TODO TODO

    #     # Get start and end indices of proteins in the system.
    #     # The assumption here is that proteins are ordered and the start residue of the next
    #     # protein is always smaller than the last residue of the previous protein.
    #     resseq = self.pdf.resSeq.to_list()
    #     p0 = resseq[0]
    #     # system_proteins contains the start and end indices of all proteins.
    #     fi_li = []
    #     fi = 0
    #     for li, p in enumerate(resseq):
    #         if p < p0:
    #             fi_li.append((fi, li-1))
    #             fi = li
    #         p0 = p
    #     fi_li.append((fi, li))

    #     self.fi_li = fi_li

    # def system_proteins(self, merge=True):
    #     """Stores information for each protein in the system using the Protein class.
    #     Returns
    #     -------
    #     proteins : list
    #         list of Protein classes. One for each different protein in the system.
    #     """
    #     if self.resolution == "martini":
    #         ref_atom = 'BB'
    #     elif self.resolution == "atomistic":
    #         ref_atom = 'CA'

    #     c = 0
    #     proteins = []
    #     # Two proteins are the same if residue number and all beads are equal between them.
    #     for values in self.fi_li:

    #         # first and last index
    #         fi = values[0]
    #         li = values[1]

    #         # Get the data for the current protein.
    #         current_protein_df = self.pdf[(self.pdf.index >= fi) & (self.pdf.index <= li)]
    #         curr_len = len(current_protein_df[current_protein_df.name == ref_atom])

    #         # curr_len = len(current_protein_df)
    #         curr_names = current_protein_df.name.to_list()

    #         new_protein = True
    #         if merge:
    #             for pc in proteins:
    #                 if curr_names == pc.beads:
    #                     pc.counter += 1
    #                     pc.dataframe.append(current_protein_df.copy())
    #                     new_protein = False

    #         if new_protein:
    #             protein = Protein(f'Protein{c}')
    #             protein.dataframe.append(current_protein_df.copy())
    #             protein.beads = curr_names
    #             protein.n_residues = curr_len
    #             protein.residues = current_protein_df.resSeq.unique()

    #             resnames = current_protein_df[current_protein_df.name == ref_atom].resName.to_numpy()
    #             protein.resnames = resnames

    #             protein.first_residue = current_protein_df.resSeq.to_list()[0]
    #             protein.last_residue = current_protein_df.resSeq.to_list()[-1]

    #             protein.resolution = self.resolution

    #             proteins.append(protein)

    #         c += 1

    #     return proteins

    # def merge_chains(self, proteins):
    #     """If the input structure file contained multiple and similar chains, ProLint
    #     will count them as separate proteins by default. This function allows you to
    #     change that. You call it on the list of proteins that you want to change and it will
    #     merge those proteins into multiple chains.
    #     Parameters
    #     ----------
    #     proteins : ProLint.Proteins
    #         Proteins object.
    #     """
    #     merged = []
    #     for protein in proteins:
    #         protein.beads = protein.beads * protein.counter
    #         protein.dataframe = [pd.concat(protein.dataframe)]
    #         # protein.n_residues = protein.n_residues #* protein.counter => n_residues remains the same
    #         protein.chains = protein.counter
    #         protein.counter = 1
    #         merged.append(protein)
    #         print ("Merged copies as chains in {}".format(protein.name))

    #     return merged

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
