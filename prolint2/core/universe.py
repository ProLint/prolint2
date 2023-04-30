import warnings

import MDAnalysis as mda

from prolint2.core.base import MacrosClass
from prolint2.core.groups import ExtendedAtomGroup
from prolint2.metrics.registries import MetricRegistry
from prolint2.core.contact_provider import ContactsProvider

warnings.filterwarnings('ignore')
class Universe(mda.Universe):
    """A subclass of MDAnalysis.Universe that adds a query and database attribute, and other useful methods."""
    def __init__(self, *args, universe=None, query=None, database=None, **kwargs):
        if universe is not None:
            if isinstance(universe, mda.Universe):
                topology = universe.filename
                trajectory = universe.trajectory.filename
                super().__init__(topology, trajectory)
            else:
                raise TypeError("universe argument should be an instance of mda.Universe")
        else:
            super().__init__(*args, **kwargs)

        self._query = self._handle_query(query)
        self._database = self._handle_database(database)

        self.registry = MetricRegistry()
        self.contacts = ContactsProvider(self.query, self.database)

        self._add_macros()

    def _add_macros(self):
        macros_attr = MacrosClass(self)
        self.atoms.universe.add_TopologyAttr(macros_attr)
        macros_attr.set_macros_values(self.query)

    def _handle_query(self, query):
        if query is None:
            query_selection_string = "protein"
            query = self.select_atoms(query_selection_string)
        return ExtendedAtomGroup(query)

    def _handle_database(self, database):
        if database is None:
            database_selection_string = "not protein"
            database = self.select_atoms(database_selection_string)
        return ExtendedAtomGroup(database)

    @property
    def query(self):
        """The query AtomGroup. This is the group of atoms that are used as reference during contact calculation."""
        return ExtendedAtomGroup(self._query)

    @query.setter
    def query(self, new_query):
        if not isinstance(new_query, mda.AtomGroup):
            raise TypeError("query attribute must be an instance of mda.AtomGroup")
        self._query = new_query

    def update_query(self, new_query):
        """Update the query AtomGroup with a new AtomGroup."""
        self.query = new_query

    @property
    def database(self):
        """The database AtomGroup. This is the group of atoms that are used as target during contact calculation."""
        return ExtendedAtomGroup(self._database)

    @database.setter
    def database(self, new_database):
        if not isinstance(new_database, mda.AtomGroup):
            raise TypeError("database attribute must be an instance of mda.AtomGroup")
        self._database = new_database

    def update_database(self, new_database):
        """Update the database AtomGroup with a new AtomGroup."""
        self.database = new_database

    def compute_contacts(self, *args, **kwargs):
        """Compute contacts between the query and database AtomGroups."""
        return self.contacts.compute(*args, **kwargs)

    def __str__(self) -> str:
        return f"<ProLint Wrapper for {super().__str__()}>"
    
    def __repr__(self) -> str:
        return f"<ProLint Wrapper for {super().__repr__()}>"
