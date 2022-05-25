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
from .contacts import Contacts, ContactsPar


class MacrosClass(ResidueStringAttr):
    attrname = 'macros'
    singular = 'macro'

    @staticmethod
    def _gen_initial_values(n_atoms, n_residues, n_segments):
        return np.array(['other'] * n_residues, dtype=object)


class ContactRunner(object):

    def __init__(self):
        self.backend = None
        self.n_jobs = -1
        # TODO
        # add funcionalities to run analysis on HPC machines


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
        self.query = self.atoms.select_atoms('')  # returns empty atomgroup
        self.database = self.atoms.select_atoms('')  # returns empty atomgroup
        self.runner = ContactRunner()
        self.cutoff = 7
        self.contacts = None

    def get_AG(self, selection, add_filter):
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
            selection = self.atoms.select_atoms(f"index {selection.index}")
        elif isinstance(selection, np.ndarray):
            selection = self.atoms[selection]
        elif isinstance(selection, str):
            selection = self.atoms.select_atoms(selection)
        return selection.select_atoms(add_filter)

    def select_query(self, selection='all', add_filter='all'):
        self.query = self.get_AG(selection, add_filter)

    def select_database(self, selection='all', add_filter='all'):
        self.database = self.get_AG(selection, add_filter)

    def get_contacts(self):
        assert isinstance(
            self.query,
            (mda.core.groups.AtomGroup),
        ), "the query has to be an AtomGroup"
        assert isinstance(
            self.database,
            (mda.core.groups.AtomGroup),
        ), "the database has to be an AtomGroup"
        if self.runner.backend == None or self.runner.backend not in ['serial', 'parallel']:
            raise ValueError(
                "You have to select a proper backend before running the contacts routine. \n Valid options: 'serial', 'parallel'"
            )
        if self.runner.backend == 'serial':
            temp_instance = Contacts(self.atoms.universe, self.query, self.database, self.cutoff)
            temp_instance.run(verbose=True)
        elif self.runner.backend == 'parallel':
            temp_instance = ContactsPar(self.atoms.universe, self.query, self.database, self.cutoff)
            temp_instance.run(n_jobs=self.runner.n_jobs)
        self.contacts = temp_instance.contacts
