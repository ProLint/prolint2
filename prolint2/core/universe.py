r""":mod:`prolint2.core.universe`
==========================================================
:Authors: Daniel P. Ramirez & Besian I. Sejdiu
:Year: 2022
:Copyright: MIT License
"""

import warnings
from typing import Literal, get_args

import os
import numpy as np
import MDAnalysis as mda

from prolint2.core.groups import ExtendedAtomGroup
from prolint2.metrics.registries import MetricRegistry
from prolint2.core.contact_provider import ContactsProvider

from prolint2.config.units import UnitConversionFactor

import configparser

# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "../config.ini")
config.read(config_file)
parameters_config = config["Parameters"]

warnings.filterwarnings("ignore")

TimeUnitLiteral = Literal["fs", "ps", "ns", "us", "ms", "s"]

# Build VALID_UNITS from TimeUnitLiteral
VALID_UNITS = get_args(TimeUnitLiteral)


class Universe(mda.Universe):
    """
    A subclass of MDAnalysis.Universe that adds a query and database attribute, and other useful methods.

    Parameters:
    *args: Variable positional arguments for the parent class.
    universe: mda.Universe or None, optional. An existing Universe object to use as the basis for this Universe.
    query: mda.AtomGroup or None, optional. The AtomGroup used as a query during contact calculations.
    database: mda.AtomGroup or None, optional. The AtomGroup used as a database during contact calculations.
    normalize_by: Literal['counts', 'actual_time', 'time_fraction'], optional. The normalization method for time.
    units: TimeUnitLiteral, optional. The units for time conversion.
    add_lipid_types: list, optional. Additional lipid types to include in the database.

    Attributes:
    query: ExtendedAtomGroup. The query AtomGroup used as a reference during contact calculation.
    database: ExtendedAtomGroup. The database AtomGroup used as a target during contact calculation.
    units: TimeUnitLiteral. The units for time conversion.
    normalize_by: str. The normalization method for time.
    """

    def __init__(
        self,
        *args,
        universe=None,
        query=None,
        database=None,
        normalize_by: Literal[
            "counts", "actual_time", "time_fraction"
        ] = "time_fraction",
        units: TimeUnitLiteral = "us",
        add_lipid_types: list = [],
        **kwargs,
    ):
        """
        Initialize the Universe.

        Args:
        *args: Variable positional arguments for the parent class.
        universe: mda.Universe or None, optional. An existing Universe object to use as the basis for this Universe.
        query: mda.AtomGroup or None, optional. The AtomGroup used as a query during contact calculations.
        database: mda.AtomGroup or None, optional. The AtomGroup used as a database during contact calculations.
        normalize_by: Literal['counts', 'actual_time', 'time_fraction'], optional. The normalization method for time.
        units: TimeUnitLiteral, optional. The units for time conversion.
        add_lipid_types: list, optional. Additional lipid types to include in the database.

        Returns:
        None
        """
        if universe is not None:
            if isinstance(universe, mda.Universe):
                topology = universe.filename
                trajectory = universe.trajectory.filename
                super().__init__(topology, trajectory)
            else:
                raise TypeError(
                    "universe argument should be an instance of mda.Universe"
                )
        else:
            super().__init__(*args, **kwargs)

        self._query = self._handle_query(query)
        # adding additional lipid types to the database
        if add_lipid_types:
            unique_lipids = (
                parameters_config["lipid_types"] + ", " + ", ".join(add_lipid_types)
            )
            unique_lipids = np.unique(unique_lipids.split(", "))
            config.set("Parameters", "lipid_types", ", ".join(unique_lipids))
            with open(config_file, "w") as configfile:
                config.write(configfile, space_around_delimiters=True)
        self._database = self._handle_database(database)

        self.params = {
            "units": units,
            "normalizer": normalize_by,
            "unit_conversion_factor": self._handle_units(units),
            "norm_factor": self._handle_normalizer(normalize_by, units),
        }

        self.registry = MetricRegistry()

    def _handle_query(self, query):
        """
        Handle the query AtomGroup.

        Args:
        query: mda.AtomGroup or None. The AtomGroup used as a query during contact calculations.

        Returns:
        ExtendedAtomGroup. The query AtomGroup.
        """
        if query is None:
            query_selection_string = "protein"
            query = self.select_atoms(query_selection_string)
        return ExtendedAtomGroup(query)

    def _handle_database(self, database):
        """
        Handle the database AtomGroup.

        Args:
        database: mda.AtomGroup or None. The AtomGroup used as a database during contact calculations.

        Returns:
        ExtendedAtomGroup. The database AtomGroup.
        """
        if database is None:
            # defining lipid types to be included in the database
            lipid_types = parameters_config["lipid_types"].split(", ")
            not_protein_restypes = np.unique(
                self.atoms.select_atoms("not protein").residues.resnames
            )
            membrane_restypes = []
            for type in lipid_types:
                if type in not_protein_restypes:
                    membrane_restypes.append("resname " + type)
            if len(membrane_restypes) == 1:
                database_selection_string = membrane_restypes[0]
            elif len(membrane_restypes) > 1:
                database_selection_string = membrane_restypes[0]
                for type in membrane_restypes[1:]:
                    database_selection_string = (
                        database_selection_string + " or " + type
                    )
            else:
                print("There are not lipid residues in your system")
            database = self.select_atoms(database_selection_string)
        return ExtendedAtomGroup(database)

    def _handle_units(self, units):
        """
        Handle time units.

        Args:
        units: TimeUnitLiteral or str. The units for time conversion.

        Returns:
        float. The unit conversion factor.
        """
        if isinstance(units, str):
            if units in UnitConversionFactor.__members__:
                units = UnitConversionFactor[units]
            else:
                raise ValueError(
                    f"units argument must be one of {UnitConversionFactor.__members__}"
                )
        time_unit = self._set_default_time_unit()
        return UnitConversionFactor[time_unit].value / units.value

    def _handle_normalizer(self, normalize_by, units):
        """
        Handle the time normalizer.

        Args:
        normalize_by: str. The normalization method for time.
        units: TimeUnitLiteral. The units for time conversion.

        Returns:
        float. The normalization factor.
        """
        if normalize_by not in ["counts", "actual_time", "time_fraction"]:
            raise ValueError(
                "normalize_by argument must be one of ['counts', 'actual_time', 'time_fraction']"
            )
        norm_factors = {
            "counts": 1.0,
            "actual_time": float(self.trajectory.dt * self._handle_units(units)),
            "time_fraction": float(self.trajectory.dt / self.trajectory.totaltime),
        }
        return norm_factors[normalize_by]

    def _set_default_time_unit(self):
        """
        Set the default time unit.

        Returns:
        str. The default time unit.
        """
        traj_time_unit = self.trajectory.units.get("time", None)
        if traj_time_unit is None:
            warnings.warn("Trajectory time unit is not set. Assuming 'ps'.")

        return traj_time_unit if traj_time_unit is not None else "ps"

    @property
    def query(self):
        """
        The query AtomGroup.

        Returns:
        ExtendedAtomGroup. The query AtomGroup.
        """
        return ExtendedAtomGroup(self._query)

    @query.setter
    def query(self, new_query):
        """
        Set the query AtomGroup.

        Args:
        new_query: mda.AtomGroup. The new query AtomGroup.

        Returns:
        None
        """
        if not isinstance(new_query, mda.AtomGroup):
            raise TypeError("query attribute must be an instance of mda.AtomGroup")
        self._query = new_query

    def update_query(self, new_query):
        """
        Update the query AtomGroup with a new AtomGroup.

        Args:
        new_query: mda.AtomGroup. The new query AtomGroup.

        Returns:
        None
        """
        self.query = new_query

    def write_lipid_occupancies_to_bfactor(self, occupancies_dict=None, lipid_type=None):
        """
        Write lipid occupancies to `bfactor` topology attribute.

        Args:
        occupancies_dict: dict, optional. A dictionary of occupancies for each lipid type.
        lipid_name: str, optional. The name of the lipid to search for occupancies.

        Returns:
        None
        """
        if occupancies_dict is None:
            # Raise an error if no metrics have been provided
            raise ValueError("No dictionary of occupancies have been provided.")
        elif lipid_type is None:
            # Raise an error if no lipid type has been provided
            raise ValueError("No lipid type has been provided.")
        else:
            occupancy_values = []
            for res in self.residues:
                if res.resid in occupancies_dict.keys():
                    occupancy_values.extend([occupancies_dict[res.resid][lipid_type]] * res.atoms.n_atoms)
                else:
                    occupancy_values.extend([0] * res.atoms.n_atoms)
        assert len(occupancy_values) == self.atoms.n_atoms
        self.add_TopologyAttr("bfactor", occupancy_values)

    def write_metrics_to_bfactor(self, metrics_dict=None, lipid_type=None, metric_name=None):
        """
        Write metrics to `bfactor` topology attribute.

        Args:
        metrics_dict: dict, optional. A dictionary of metrics for each lipid type.
        lipid_name: str, optional. The name of the lipid to search for metrics.
        metric_name: str, optional. The name of the metric to write.

        Returns:
        None
        """
        if metrics_dict is None:
            # Raise an error if no metrics have been provided
            raise ValueError("No dictionary of metrics have been provided.")
        elif lipid_type is None:
            # Raise an error if no lipid type has been provided
            raise ValueError("No lipid type has been provided.")
        elif metric_name is None:
            # Raise an error if no metric name has been provided
            raise ValueError("No metric name has been provided.")
        else:
            metric_values = []
            for res in self.residues:
                if res.resid in metrics_dict.keys() and metric_name in metrics_dict[res.resid][lipid_type].keys():
                    metric_values.extend([metrics_dict[res.resid][lipid_type][metric_name]] * res.atoms.n_atoms)
                else:
                    metric_values.extend([0] * res.atoms.n_atoms)
        assert len(metric_values) == self.atoms.n_atoms
        # normalize values to be between 0 and 1
        metric_values = (np.array(metric_values) - np.min(metric_values)) * 100 / (np.max(metric_values) - np.min(metric_values))
        self.add_TopologyAttr("bfactor", metric_values)



    @property
    def database(self):
        """
        The database AtomGroup.

        Returns:
        ExtendedAtomGroup. The database AtomGroup.
        """
        return ExtendedAtomGroup(self._database)

    @database.setter
    def database(self, new_database):
        """
        Set the database AtomGroup.

        Args:
        new_database: mda.AtomGroup. The new database AtomGroup.

        Returns:
        None
        """
        if not isinstance(new_database, mda.AtomGroup):
            raise TypeError("database attribute must be an instance of mda.AtomGroup")
        self._database = new_database

    def update_database(self, new_database):
        """
        Update the database AtomGroup with a new AtomGroup.

        Args:
        new_database: mda.AtomGroup. The new database AtomGroup.

        Returns:
        None
        """
        self.database = new_database

    def compute_contacts(self, *args, **kwargs):
        """
        Compute contacts between the query and database AtomGroups.

        Args:
        *args: Additional positional arguments for ContactsProvider.compute.
        **kwargs: Additional keyword arguments for ContactsProvider.compute.

        Returns:
        result: The result of the contact computation.
        """
        contacts_provider = ContactsProvider(
            self.query, self.database, params=self.params
        )
        return contacts_provider.compute(*args, **kwargs)

    def load_contacts_from_file(self, file, *args, **kwargs):
        """
        Load contacts from a file.

        Args:
        file: str. The path to the file containing contacts.
        *args: Additional positional arguments for ContactsProvider.load_from_file.
        **kwargs: Additional keyword arguments for ContactsProvider.load_from_file.

        Returns:
        result: The result of the contact loading operation.
        """
        contacts_provider = ContactsProvider(
            self.query, self.database, params=self.params
        )
        return contacts_provider.load_from_file(file, *args, **kwargs)

    @property
    def units(self):
        """
        The units of the trajectory time.

        Returns:
        TimeUnitLiteral. The units for time conversion.
        """
        return self.params["units"]

    @units.setter
    def units(self, new_units):
        """
        Set the units for time conversion.

        Args:
        new_units: TimeUnitLiteral or str. The new units for time conversion.

        Returns:
        None
        """
        self.params["unit_conversion_factor"] = self._handle_units(new_units)
        self.params["units"] = new_units
        self.params["norm_factor"] = self._handle_normalizer(
            self.params["normalizer"], new_units
        )

    @property
    def normalize_by(self):
        """
        The normalizer of the trajectory time.

        Returns:
        str. The normalization method for time.
        """
        return self.params["normalizer"]

    @normalize_by.setter
    def normalize_by(self, new_normalizer):
        """
        Set the normalization method for time.

        Args:
        new_normalizer: str. The new normalization method for time.

        Returns:
        None
        """
        self.params["norm_factor"] = self._handle_normalizer(
            new_normalizer, self.params["units"]
        )
        self.params["normalizer"] = new_normalizer

    def __str__(self) -> str:
        """
        Return a string representation of the object.

        Returns:
        --------
        str
            A string representation of the object.
        """
        return f"<ProLint Wrapper for {super().__str__()}>"

    def __repr__(self) -> str:
        """
        Return a detailed string representation of the object.

        Returns:
        --------
        str
            A detailed string representation of the object.
        """
        return f"<ProLint Wrapper for {super().__repr__()}>"
