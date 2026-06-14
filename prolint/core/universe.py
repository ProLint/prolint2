"""ProLint Universe module.

This module provides the main entry point for ProLint, extending MDAnalysis
Universe with biomolecular interaction analysis capabilities.
"""

import logging
import warnings
from typing import Literal, Optional, Any, Union

import MDAnalysis as mda

from prolint.core.groups import ExtendedAtomGroup
from prolint.core.contact_provider import ContactsProvider
from prolint.core.replica_detection import detect_replicas, get_replica_atoms

from prolint.config.units import UnitConversionFactor
from prolint.config.enums import TimeUnit, NormalizationMethod

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)

TimeUnitLiteral = Literal["fs", "ps", "ns", "us", "ms", "s"]
NormalizationLiteral = Literal["counts", "actual_time"]

VALID_UNITS = tuple(unit.value for unit in TimeUnit)


class Universe(mda.Universe):
    """ProLint Universe for biomolecular interaction analysis.

    Extends MDAnalysis Universe with specialized methods for computing
    and analyzing contacts between molecular groups (e.g., protein-lipid,
    protein-ligand interactions).

    Parameters
    ----------
    *args : tuple
        Positional arguments passed to MDAnalysis Universe (topology, trajectory).
    universe : mda.Universe, optional
        Existing MDAnalysis Universe to wrap. If provided, topology and
        trajectory are extracted from this universe.
    query : mda.AtomGroup, optional
        Atoms to analyze (e.g., protein). Default: ``"protein"`` selection.
    database : mda.AtomGroup, optional
        Reference atoms for contact detection (e.g., lipids).
        Default: ``"not protein"`` selection.
    normalize_by : {"counts", "actual_time"}, default="counts"
        Normalization method for contact durations.
        - ``"counts"``: Duration in frame counts
        - ``"actual_time"``: Duration in time units
    units : {"fs", "ps", "ns", "us", "ms", "s"}, default="us"
        Time units for analysis results.
    **kwargs : dict
        Additional keyword arguments passed to MDAnalysis Universe.

    Examples
    --------
    Basic usage with topology and trajectory files:

    >>> from prolint import Universe
    >>> u = Universe("topology.gro", "trajectory.xtc")
    >>> print(f"Loaded {u.trajectory.n_frames} frames")

    From an existing MDAnalysis Universe:

    >>> import MDAnalysis as mda
    >>> mda_u = mda.Universe("topology.gro", "trajectory.xtc")
    >>> u = Universe(universe=mda_u)

    With custom selections:

    >>> u = Universe("topology.gro", "trajectory.xtc")
    >>> u.query = u.select_atoms("protein and name CA")
    >>> u.database = u.select_atoms("resname POPC POPE CHOL")

    Compute contacts:

    >>> contacts = u.compute_contacts(cutoff=7.0)
    >>> result = contacts.analyze("timeseries", database_type="CHOL")

    See Also
    --------
    ComputedContacts : Result object from compute_contacts
    ExtendedAtomGroup : Enhanced atom group with additional properties
    """

    def __init__(
        self,
        *args: Any,
        universe: Optional[mda.Universe] = None,
        query: Optional[mda.AtomGroup] = None,
        database: Optional[mda.AtomGroup] = None,
        normalize_by: Union[NormalizationMethod, NormalizationLiteral] = "counts",
        units: Union[TimeUnit, TimeUnitLiteral] = "us",
        **kwargs: Any,
    ) -> None:
        if universe is not None:
            if isinstance(universe, mda.Universe):
                topology = universe.filename
                trajectory = universe.trajectory.filename
                logger.info("Loading from existing Universe: %s", topology)
                super().__init__(topology, trajectory)
            else:
                raise TypeError(
                    "universe argument should be an instance of mda.Universe"
                )
        else:
            if args:
                logger.info("Loading topology: %s", args[0])
                if len(args) > 1:
                    logger.info("Loading trajectory: %s", args[1])
            super().__init__(*args, **kwargs)

        self._query = self._handle_query(query)
        self._database = self._handle_database(database)

        self.params = {
            "units": units,
            "normalizer": normalize_by,
            "unit_conversion_factor": self._handle_units(units),
            "norm_factor": self._handle_normalizer(normalize_by, units),
        }

        # Log universe summary
        logger.info(
            "Universe loaded: %d frames, %d atoms, %d residues",
            self.trajectory.n_frames,
            self.atoms.n_atoms,
            self.residues.n_residues,
        )
        logger.debug("Units: %s, Normalization: %s", units, normalize_by)

    def _handle_query(self, query: Optional[mda.AtomGroup]) -> ExtendedAtomGroup:
        """Initialize query atom group.

        Parameters
        ----------
        query : mda.AtomGroup, optional
            User-provided query selection.

        Returns
        -------
        ExtendedAtomGroup
            Wrapped query atom group.
        """
        if query is None:
            query_selection_string = "protein"
            query = self.select_atoms(query_selection_string)
            logger.info(
                "Query selection: '%s' (%d atoms)",
                query_selection_string,
                query.n_atoms,
            )
        else:
            logger.info("Query selection: custom (%d atoms)", query.n_atoms)
        return ExtendedAtomGroup(query)

    def _handle_database(self, database: Optional[mda.AtomGroup]) -> ExtendedAtomGroup:
        """Initialize database atom group.

        Parameters
        ----------
        database : mda.AtomGroup, optional
            User-provided database selection.

        Returns
        -------
        ExtendedAtomGroup
            Wrapped database atom group.
        """
        if database is None:
            database_selection_string = "not protein"
            database = self.select_atoms(database_selection_string)
            logger.info(
                "Database selection: '%s' (%d atoms)",
                database_selection_string,
                database.n_atoms,
            )
        else:
            logger.info("Database selection: custom (%d atoms)", database.n_atoms)
        return ExtendedAtomGroup(database)

    def _handle_units(self, units: Union[TimeUnit, TimeUnitLiteral, str]) -> float:
        """Convert time units to conversion factor.

        Parameters
        ----------
        units : TimeUnit or str
            Time unit specification.

        Returns
        -------
        float
            Conversion factor from trajectory time unit to specified unit.

        Raises
        ------
        ValueError
            If units is not a valid time unit.
        """
        # Convert TimeUnit enum to string if needed
        if isinstance(units, TimeUnit):
            units_str = units.value
        else:
            units_str = units

        units_enum: UnitConversionFactor
        if units_str in UnitConversionFactor.__members__:
            units_enum = UnitConversionFactor[units_str]
        else:
            raise ValueError(
                f"units argument must be one of {list(UnitConversionFactor.__members__.keys())}"
            )
        time_unit = self._set_default_time_unit()
        return UnitConversionFactor[time_unit].value / units_enum.value

    def _handle_normalizer(
        self,
        normalize_by: Union[NormalizationMethod, NormalizationLiteral],
        units: Union[TimeUnit, TimeUnitLiteral, str],
    ) -> float:
        """Compute normalization factor for contact durations.

        Parameters
        ----------
        normalize_by : NormalizationMethod or str
            Normalization method.
        units : TimeUnit or str
            Time units for normalization.

        Returns
        -------
        float
            Normalization factor.

        Raises
        ------
        ValueError
            If normalize_by is not a valid normalization method.
        """
        # Convert enum to string if needed
        if isinstance(normalize_by, NormalizationMethod):
            normalize_str = normalize_by.value
        else:
            normalize_str = normalize_by

        valid_methods = [m.value for m in NormalizationMethod]
        if normalize_str not in valid_methods:
            raise ValueError(f"normalize_by argument must be one of {valid_methods}")
        norm_factors = {
            "counts": 1.0,
            "actual_time": float(self.trajectory.dt * self._handle_units(units)),
        }
        return norm_factors[normalize_str]

    def _set_default_time_unit(self) -> str:
        """Get default time unit from trajectory.

        Returns
        -------
        str
            Time unit string (defaults to "ps" if not set in trajectory).
        """
        traj_time_unit = self.trajectory.units.get("time", None)
        if traj_time_unit is None:
            warnings.warn("Trajectory time unit is not set. Assuming 'ps'.")

        return traj_time_unit if traj_time_unit is not None else "ps"

    @property
    def query(self) -> ExtendedAtomGroup:
        """Query atom group for analysis.

        Returns
        -------
        ExtendedAtomGroup
            The query atoms (typically protein).
        """
        return self._query

    @query.setter
    def query(self, new_query: mda.AtomGroup) -> None:
        """Set query atom group.

        Parameters
        ----------
        new_query : mda.AtomGroup
            New query atom selection.

        Raises
        ------
        TypeError
            If new_query is not an AtomGroup.
        """
        if not isinstance(new_query, mda.AtomGroup):
            raise TypeError("query attribute must be an instance of mda.AtomGroup")
        self._query = (
            new_query
            if isinstance(new_query, ExtendedAtomGroup)
            else ExtendedAtomGroup(new_query)
        )

    def update_query(self, new_query: mda.AtomGroup) -> None:
        """Update query atom group.

        Parameters
        ----------
        new_query : mda.AtomGroup
            New query atom selection.
        """
        self.query = new_query

    @property
    def database(self) -> ExtendedAtomGroup:
        """Database atom group for contact detection.

        Returns
        -------
        ExtendedAtomGroup
            The database atoms (typically lipids or ligands).
        """
        return self._database

    @database.setter
    def database(self, new_database: mda.AtomGroup) -> None:
        """Set database atom group.

        Parameters
        ----------
        new_database : mda.AtomGroup
            New database atom selection.

        Raises
        ------
        TypeError
            If new_database is not an AtomGroup.
        """
        if not isinstance(new_database, mda.AtomGroup):
            raise TypeError("database attribute must be an instance of mda.AtomGroup")
        self._database = (
            new_database
            if isinstance(new_database, ExtendedAtomGroup)
            else ExtendedAtomGroup(new_database)
        )

    def update_database(self, new_database: mda.AtomGroup) -> None:
        """Update database atom group.

        Parameters
        ----------
        new_database : mda.AtomGroup
            New database atom selection.
        """
        self.database = new_database

    def compute_contacts(
        self,
        *args: Any,
        replica: Optional[str] = None,
        **kwargs: Any,
    ) -> Any:
        """Compute contacts between query and database atom groups.

        Detects distance-based contacts between query (e.g., protein) and
        database (e.g., lipids) atoms across the trajectory.

        Parameters
        ----------
        cutoff : float, default=7.0
            Distance cutoff in Angstroms for contact detection.
        start : int, optional
            First frame to process. Default: first frame.
        stop : int, optional
            Last frame to process (exclusive). Default: last frame.
        step : int, default=1
            Frame step size.
        replica : str, optional
            Replica ID to analyze (e.g., 'A', 'B'). Required when multiple
            replicas with repeated residue IDs are detected. Use
            ``detect_replicas()`` to see available replicas.
        **kwargs : dict
            Additional keyword arguments passed to ContactsProvider.

        Returns
        -------
        ComputedContacts
            Object containing contact data and analysis methods.

        Raises
        ------
        ValueError
            If multiple replicas with repeated residue IDs are detected and
            ``replica`` is not specified.

        Examples
        --------
        Compute contacts for a single-replica system:

        >>> contacts = universe.compute_contacts(cutoff=7.0)

        For multi-replica systems, first detect replicas:

        >>> from prolint.core.replica_detection import detect_replicas
        >>> result = detect_replicas(universe.query)
        >>> print(f"Found {result.n_replicas} replicas")
        >>> for info in result.replica_info:
        ...     print(f"  Replica {info.replica_id}: {info.n_residues} residues")

        Then select a specific replica for analysis:

        >>> contacts = universe.compute_contacts(cutoff=7.0, replica='A')

        See Also
        --------
        ComputedContacts : Result object with analysis methods
        detect_replicas : Detect replicas in query selection
        """
        if "cutoff" not in kwargs:
            warnings.warn("No cutoff provided. Using default cutoff of 7.0 Angstroms.")
            kwargs["cutoff"] = 7.0

        # Detect replicas in query selection
        replica_result = detect_replicas(self.query)

        if replica_result.n_replicas > 1:
            if replica_result.has_repeated_resids:
                # Multiple replicas with repeated residue IDs - must select one
                if replica is None:
                    replica_ids = [info.replica_id for info in replica_result.replica_info]
                    raise ValueError(
                        f"Detected {replica_result.n_replicas} replicas with repeated "
                        f"residue IDs. Please specify which replica to analyze using "
                        f"replica='{replica_ids[0]}' (or one of: {replica_ids})."
                    )

                # Update query to selected replica
                self.query = get_replica_atoms(replica_result, replica)
                replica_info = next(
                    info for info in replica_result.replica_info if info.replica_id == replica
                )
                logger.info(
                    f"Analyzing replica {replica} ({replica_info.n_residues} residues)"
                )
            elif replica is not None:
                # Different sequences but user selected a specific replica
                self.query = get_replica_atoms(replica_result, replica)
                replica_info = next(
                    info for info in replica_result.replica_info if info.replica_id == replica
                )
                logger.info(
                    f"Analyzing replica {replica} ({replica_info.n_residues} residues)"
                )
            else:
                logger.info(
                    f"Analyzing {replica_result.n_replicas} replicas with unique residue IDs"
                )

        contacts_provider = ContactsProvider(
            self.query,
            self.database,
            params=self.params,
        )
        return contacts_provider.compute(*args, **kwargs)

    @property
    def units(self) -> Union[TimeUnit, str]:
        """Current time units for analysis results.

        Returns
        -------
        TimeUnit or str
            Current time unit setting.
        """
        return self.params["units"]

    @units.setter
    def units(self, new_units: Union[TimeUnit, TimeUnitLiteral, str]) -> None:
        """Set time units for analysis results.

        Parameters
        ----------
        new_units : TimeUnit or str
            New time unit (fs, ps, ns, us, ms, s).
        """
        self.params["unit_conversion_factor"] = self._handle_units(new_units)
        self.params["units"] = new_units
        self.params["norm_factor"] = self._handle_normalizer(
            self.params["normalizer"], new_units
        )

    @property
    def normalize_by(self) -> Union[NormalizationMethod, str]:
        """Current normalization method for contact durations.

        Returns
        -------
        NormalizationMethod or str
            Current normalization setting.
        """
        return self.params["normalizer"]

    @normalize_by.setter
    def normalize_by(
        self, new_normalizer: Union[NormalizationMethod, NormalizationLiteral]
    ) -> None:
        """Set normalization method for contact durations.

        Parameters
        ----------
        new_normalizer : NormalizationMethod or str
            New normalization method (counts or actual_time).
        """
        self.params["norm_factor"] = self._handle_normalizer(
            new_normalizer, self.params["units"]
        )
        self.params["normalizer"] = new_normalizer

    def __str__(self) -> str:

        return f"<ProLint Wrapper for {super().__str__()}>"

    def __repr__(self) -> str:
        return f"<ProLint Wrapper for {super().__repr__()}>"
