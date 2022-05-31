r"""UFCC base classes --- :mod:`ufcc.ufcc`
======================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import ResidueStringAttr
from .contacts import Contacts


class MacrosClass(ResidueStringAttr):
    """
    Class to add the *macros* metadata. 
    
    The *macros* metadata is an aditional label to each reidue in the system,
    that is going to be useful for the selection of the query and the database groups.
  
    If the residue is included in the :class:`MembraneDatabase` group, then the *macro*
    metadata will be setted as `membrane`; if the residue is included in the
    :class:`QueryProteins` group then the *macro* metadata will be settep with 
    **protein** or **protein#** depending on the number of segments (or chains) in the system; 
    otherwise the *macro* metadata will be setted as **other**.

    .. warning::

        The identification of the different proteins in the system will be done using one of two
        approaches:

        i. If the format files used include segment (or chain) information, then the `macro` metadata will
        be setted with the name specified in each segment (or chain).
        ii. If the format files used do not include this information (i.e. *gro* format file) then :class:`UFCC`
        will assume that proteins are ordered and the start residue of the next protein is always smaller than 
        the last residue of the previous protein. 
   
    Example
    -------
    All these assignnation are done automatically by the :class:`UFCC`, so you do not need to use this 
    class for anything. But you can access the information of the *macros metadata* as follows::

        from ufcc import UFCC
        target_system = UFCC('coordinates.gro', 'trajectory.xtc')

        target_system.query.selected.residues.macros

    And you will get an uni-dimensional numpy array with same amount of values as residues selected in the **query**
    group and the **macro** of each residue. You can do the same for your **database** group.

    """
    attrname = 'macros'
    singular = 'macro'

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments):
        return np.array(['other'] * n_residues, dtype=object)


class UFCC(object):
    """Base class for managing the distance-based contacts calculation routines done by the **UFCC** 
    package. It reads a structure/topology file and a trajectory file in any of the MDAnalysis-supported
    formats.

    Parameters
    ----------
    structure: any MDAnalysis-supported structure/topology file.

    trajectory : any MDAnalysis-supported structure/topology file.

    add_lipid_types : list
        list of strings with the residue name of lipids not included in the **UFCC** list of supported residues.

    Attributes
    ----------
    atoms : AtomGroup
        MDAnalysis AtomGroup with all the atoms in the system.
    residues : ResidueGroup
        MDAnalysis ResidueGroup with all the residues in the system.
    list_of_macros : list
        All the availables macros to use during the selection of the query/database groups.
    query : :class:`QueryProteins`
        **Query** group to use during the contacts calculation.
    database : :class:`MembraneDatabase`
        **Database** group to use during the contacts calculation.
    contacts : :class:`.contacts.Contacts`
        Contacts object to run and analyze the distance-based contacts results.
    """

    def __init__(self, structure, trajectory, add_lipid_types=[]):
        # wrapping some basic MDAnalysis groups
        self.atoms = mda.Universe(structure, trajectory).atoms
        self.residues = self.atoms.residues
        self.atoms.universe.add_TopologyAttr('macros')

        # adding the macros to the membrane residues
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

        # adding the macros to the protein residues
        protein_sel = 'protein'
        # First possibility: we can access segment(chain) information from the Universe.
        if len(self.atoms.select_atoms(protein_sel).segments) > 1 and self.atoms.select_atoms(
                protein_sel).segments.n_atoms == self.atoms.select_atoms(protein_sel).n_atoms:
            for segment_idx in range(len(self.atoms.select_atoms(protein_sel).segments)):
                self.atoms.select_atoms(
                    protein_sel).segments[segment_idx].residues.macros = 'Protein' + str(segment_idx)
        # Second possibility: the assumption here is that proteins are ordered and the start residue of the next
        # protein is always smaller than the last residue of the previous protein.
        else:
            # Get start and end indices of proteins in the system.
            resseq = self.atoms.select_atoms(protein_sel).residues.resindices
            p0 = resseq[0]
            # first and last index
            fi_li = []
            fi = 0
            for li, p in enumerate(resseq):
                if p < p0:
                    fi_li.append((fi, li - 1))
                    fi = li
                p0 = p
            fi_li.append((fi, li))

            for idx, values in enumerate(fi_li):
                fi = values[0]
                li = values[1]
                self.atoms.select_atoms(protein_sel).residues[list(range(fi, li +
                                                                         1))].residues.macros = 'Protein' + str(idx)

        # TODO
        # Add merge chains and options to change the name of the proteins.

        self.atoms.select_atoms(membrane_sel).residues.macros = 'membrane'
        self.list_of_macros = list(np.unique(self.atoms.residues.macros))
        self.query = QueryProteins(self.atoms.select_atoms(protein_sel))
        self.database = MembraneDatabase(self.atoms.select_atoms(membrane_sel))
        self.contacts = Contacts(self.query, self.database)

    def __str__(self):
        return "Base class to handle the contacts calculation in UFCC."

    def __repr__(self):
        return "Base class to handle the contacts calculation in UFCC."


class BasicGroup(object):
    """ 
    Basic class to be heritaged for the :class:`MembraneDatabase` and :class:`QueryProteins`
    classes in order to handle the **database** and **query** groups respectivaly.

    Attributes
    ----------
    selected : AtomGroup
        An MDAnalysis AtomGroup object that includes the atoms 
        that will be used as database/query for the contacts calculation.
    whole : AtomGroup
        An MDAnalysis AtomGroup object including all the atoms from where the selections
        can be done to define the database/query atoms for the contacts calculation. 
        The **selected** attribute will be always a subset of the **whole**. 
    """

    def __init__(self, whole):
        self.selected = whole
        self.whole = whole

    def select(self, selection='all'):
        """
        Cast an MDAnalysis.Atom, MDAnalysis.Residue, MDAnalysis.ResidueGroup, or str syntax 
        from the **whole** AtomGroup to the **selected** AtomGroup.

        Parameters
        ----------
        selection: MDAnalysis.Atom, MDAnalysis.Residue,  MDAnalysis.ResidueGroup or str
            atoms to cast
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
        self.selected = selection


class MembraneDatabase(BasicGroup):
    """
    Class to handle the membrane **database** group. 
    
    It heritages all atributes and methods from the 
    :class:`BasicGroup` class, and includes some new ones that
    are specific for the membrane **database** group.
    """

    def __init__(self, whole):
        super().__init__(whole)

    def lipid_types(self):
        """Get the names of all the lipids that will be analyzed.

        Returns
        -------
        array of lipid names
        """
        if not isinstance(self.selected, mda.core.groups.AtomGroup):
            return np.array([])
        else:
            return np.unique(self.selected.resnames)

    def lipid_count(self):
        """Get the name and count of each lipid that will be analyzed.

        Returns
        -------
        dictionary
            key:value corresponds to lipid_name:count.
        """
        lc = {}
        lipids = self.lipid_types()
        for lipid in lipids:
            lc[lipid] = len(self.selected.residues[self.selected.residues.resnames == lipid])
        return lc

    def __str__(self):
        if not isinstance(self.selected, mda.core.groups.AtomGroup):
            return "<ufcc.MembraneDatabase containing 0 atoms>"
        else:
            return "<ufcc.MembraneDatabase containing {} atoms>".format(self.selected.n_atoms)

    def __repr__(self):
        if not isinstance(self.selected, mda.core.groups.AtomGroup):
            return "<ufcc.MembraneDatabase containing 0 atoms>"
        else:
            return "<ufcc.MembraneDatabase containing {} atoms>".format(self.selected.n_atoms)


class QueryProteins(BasicGroup):
    """
    Class to handle the **query** proteins group. 
    
    It heritages all atributes and methods from the 
    :class:`BasicGroup` class, and includes some new ones that
    are specific for the **query** proteins group.
    """

    def __init__(self, whole):
        super().__init__(whole)

    def list_proteins(self):
        """Get the labels of all the proteins that will be analyzed.

        Returns
        -------
        array of protein labels
        """
        return np.unique(self.whole.residues.macros)

    def __str__(self):
        if not isinstance(self.selected, mda.core.groups.AtomGroup):
            return "<ufcc.QueryProteins containing 0 atoms>"
        else:
            return "<ufcc.QueryProteins containing {} atoms>".format(self.selected.n_atoms)

    def __repr__(self):
        if not isinstance(self.selected, mda.core.groups.AtomGroup):
            return "<ufcc.QueryProteins containing 0 atoms>"
        else:
            return "<ufcc.QueryProteins containing {} atoms>".format(self.selected.n_atoms)
