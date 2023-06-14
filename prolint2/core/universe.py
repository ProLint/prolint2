import warnings
from typing import Literal, get_args

import MDAnalysis as mda

from prolint2.core.base import MacrosClass
from prolint2.core.groups import ExtendedAtomGroup
from prolint2.metrics.registries import MetricRegistry
from prolint2.core.contact_provider import ContactsProvider

from prolint2.config.units import UnitConversionFactor

warnings.filterwarnings('ignore')

TimeUnitLiteral = Literal['fs', 'ps', 'ns', 'us', 'ms', 's']

# Build VALID_UNITS from TimeUnitLiteral
VALID_UNITS = get_args(TimeUnitLiteral)

class Universe(mda.Universe):
    """A subclass of MDAnalysis.Universe that adds a query and database attribute, and other useful methods."""
    def __init__(self, *args, universe=None, query=None, database=None, normalize_by: Literal['counts', 'actual_time', 'time_fraction'] = 'time_fraction', units: TimeUnitLiteral = 'us', **kwargs):
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

        self.params = {
            'units': units,
            'normalizer': normalize_by,
            'unit_conversion_factor': self._handle_units(units),
            'norm_factor': self._handle_normalizer(normalize_by, units)
        }

        self.registry = MetricRegistry()

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

    def _handle_units(self, units):
        if isinstance(units, str):
            if units in UnitConversionFactor.__members__:
                units = UnitConversionFactor[units]
            else:
                raise ValueError(f"units argument must be one of {UnitConversionFactor.__members__}")
        time_unit = self._set_default_time_unit()
        return UnitConversionFactor[time_unit].value / units.value

    def _handle_normalizer(self, normalize_by, units):
        if normalize_by not in ['counts', 'actual_time', 'time_fraction']:
            raise ValueError("normalize_by argument must be one of ['counts', 'actual_time', 'time_fraction']")
        norm_factors = {
            'counts': 1.0,
            'actual_time': float(self.trajectory.dt * self._handle_units(units)),
            'time_fraction': float(self.trajectory.dt / self.trajectory.totaltime)
        }
        return norm_factors[normalize_by]

    def _set_default_time_unit(self):
        traj_time_unit = self.trajectory.units.get('time', None)
        if traj_time_unit is None:
            warnings.warn("Trajectory time unit is not set. Assuming 'ps'.")

        return traj_time_unit if traj_time_unit is not None else 'ps'

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
        contacts_provider = ContactsProvider(self.query, self.database, params=self.params)
        return contacts_provider.compute(*args, **kwargs)

    @property
    def units(self):
        """The units of the trajectory time."""
        return self.params['units']

    @units.setter
    def units(self, new_units):
        self.params['unit_conversion_factor'] = self._handle_units(new_units)
        self.params['units'] = new_units
        self.params['norm_factor'] = self._handle_normalizer(self.params['normalizer'], new_units)
        
    @property
    def normalize_by(self):
        """The normalizer of the trajectory time."""
        return self.params['normalizer']

    @normalize_by.setter
    def normalize_by(self, new_normalizer):
        self.params['norm_factor'] = self._handle_normalizer(new_normalizer, self.params['units'])
        self.params['normalizer'] = new_normalizer

    def __str__(self) -> str:
        return f"<ProLint Wrapper for {super().__str__()}>"
    
    def __repr__(self) -> str:
        return f"<ProLint Wrapper for {super().__repr__()}>"
