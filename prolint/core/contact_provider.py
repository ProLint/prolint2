"""Contact provider module.

This module provides classes for computing and managing contact data
between molecular groups.
"""

import logging
from collections import defaultdict
from typing import Callable, Literal, TYPE_CHECKING

from prolint.computers.contacts import SerialContacts

from prolint.contacts.exact_contacts import ExactContacts

from prolint.contacts.base import BaseContactStore

from prolint.config.units import DEFAULT_SIM_PARAMS

if TYPE_CHECKING:
    from prolint.analysis.base import AnalysisResult

logger = logging.getLogger(__name__)


class ComputedContacts:
    """Container for computed contact data with analysis methods.

    This class wraps contact computation results and provides methods
    for analyzing contacts, computing metrics, and performing set operations.

    Parameters
    ----------
    contact_strategy_instance : BaseContactStore
        Contact storage instance containing computed contacts.
    provider : ContactsProvider
        The provider that created this instance.

    Examples
    --------
    >>> contacts = universe.compute_contacts(cutoff=7.0)
    >>> result = contacts.analyze("timeseries", database_type="CHOL")

    Compute metrics:

    >>> occupancy = contacts.compute_metric("occupancy", target_resname="CHOL")

    Set operations:

    >>> common = contacts1 + contacts2  # Intersection
    >>> unique = contacts1 - contacts2  # Difference

    See Also
    --------
    Universe.compute_contacts : Method that creates ComputedContacts
    AnalysisRegistry : Registry of available analysis types
    """

    def __init__(
        self, contact_strategy_instance: BaseContactStore, provider: "ContactsProvider"
    ):
        self._contact_strategy = contact_strategy_instance
        self.provider = provider

    def compute_metric(self, metric: str, target_resname=None):
        """Compute a metric for contacts.

        Parameters
        ----------
        metric : str
            Metric to compute. Options: "occupancy", "mean", "max", "sum".
        target_resname : str, optional
            Filter by residue name (e.g., "CHOL", "POPC").

        Returns
        -------
        dict
            Metric values organized by residue ID.

        Examples
        --------
        >>> occupancy = contacts.compute_metric("occupancy", target_resname="CHOL")
        >>> mean_duration = contacts.compute_metric("mean")
        """

        return self._contact_strategy.compute(metric, target_resname=target_resname)

    def apply_function(self, func: Callable, target_resname=None):
        """Apply a custom function to contact data.

        Parameters
        ----------
        func : callable
            Function to apply to contact durations.
        target_resname : str, optional
            Filter by residue name.

        Returns
        -------
        dict
            Function results organized by residue ID.
        """
        return self._contact_strategy.apply_function(
            func, target_resname=target_resname
        )

    @property
    def contacts(self):
        """Raw contact data.

        Returns
        -------
        dict
            Contact data organized by residue and database molecule.
        """
        return self._contact_strategy.contacts

    @property
    def contact_frames(self):
        """Frame indices where contacts occur.

        Returns
        -------
        dict
            Nested dict mapping residue_id -> database_id -> list of frame indices.
        """
        return self._contact_strategy.contact_frames

    @property
    def norm_factor(self) -> float:
        """Normalization factor for duration calculations.

        Returns
        -------
        float
            Factor used to normalize contact durations.
        """
        return self._contact_strategy.norm_factor

    def intersection(self, other: "ComputedContacts") -> "ComputedContacts":
        """Find contacts common to both ComputedContacts objects.

        Parameters
        ----------
        other : ComputedContacts
            Another computed contacts object.

        Returns
        -------
        ComputedContacts
            New object containing only contacts present in both.

        Examples
        --------
        >>> common = contacts1.intersection(contacts2)
        >>> # Or using operator
        >>> common = contacts1 + contacts2
        """
        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, database_ids in self.contact_frames.items():
            for database_id in database_ids:
                if database_id in other.contact_frames[residue_id]:
                    result_data[residue_id][database_id] = other.contact_frames[
                        residue_id
                    ][database_id]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(
            self.provider.query.universe, result_data
        )
        contact_instances.norm_factor = self.provider.params.get("norm_factor", 1)
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def difference(self, other: "ComputedContacts") -> "ComputedContacts":
        """Find contacts in self but not in other.

        Parameters
        ----------
        other : ComputedContacts
            Another computed contacts object.

        Returns
        -------
        ComputedContacts
            New object containing contacts unique to self.

        Examples
        --------
        >>> unique = contacts1.difference(contacts2)
        >>> # Or using operator
        >>> unique = contacts1 - contacts2
        """
        result_data = defaultdict(lambda: defaultdict(list))

        for residue_id, database_ids in self.contact_frames.items():
            for database_id in database_ids:
                if database_id not in other.contact_frames[residue_id]:
                    result_data[residue_id][database_id] = self.contact_frames[
                        residue_id
                    ][database_id]

        # Create a new instance of the contact strategy class
        contact_instances = self._contact_strategy.__class__(
            self.provider.query.universe, result_data
        )
        contact_instances.run()

        return ComputedContacts(contact_instances, self.provider)

    def __add__(self, other: "ComputedContacts") -> "ComputedContacts":
        """Intersection operator. See :meth:`intersection`."""
        return self.intersection(other)

    def __sub__(self, other: "ComputedContacts") -> "ComputedContacts":
        """Difference operator. See :meth:`difference`."""
        return self.difference(other)

    def analyze(self, analysis_type: str, **kwargs) -> "AnalysisResult":
        """Run an analysis on computed contacts.

        Parameters
        ----------
        analysis_type : str
            Type of analysis to run. Options:
            - "timeseries": Contact counts over time
            - "database_contacts": Per-molecule contact timeline
            - "kinetics": Binding kinetics and residence times
            - "density_map": 2D spatial density
            - "radial_density": Radial density profile
            - "shared_contacts": Residue contact correlations
            - "distances": Distance distributions
            - "atom_distances": Atom-level distances
            - "metrics": Per-residue metrics
        **kwargs : dict
            Analysis-specific parameters.

        Returns
        -------
        AnalysisResult
            Result object containing analysis data.

        Examples
        --------
        >>> result = contacts.analyze("timeseries", database_type="CHOL")
        >>> result = contacts.analyze("kinetics", query_residue=42, mode="accumulated")
        >>> result = contacts.analyze("metrics", metric="occupancy")

        See Also
        --------
        AnalysisRegistry : Registry of available analysis types
        """
        from prolint.analysis import AnalysisRegistry

        analysis = AnalysisRegistry.create(
            analysis_type,
            self.provider.query.universe,
            self,
        )
        return analysis.run(**kwargs)


class ContactsProvider:
    """Orchestrates contact computation between atom groups.

    This class manages the contact detection process, coordinating
    between different computation strategies and contact storage methods.

    Parameters
    ----------
    query : ExtendedAtomGroup
        Query atoms (e.g., protein).
    database : ExtendedAtomGroup
        Database atoms (e.g., lipids).
    params : dict, optional
        Computation parameters including units and normalization.
    compute_strategy : {"default"}, default="default"
        Contact computation strategy to use.

    See Also
    --------
    Universe.compute_contacts : High-level interface using this provider
    ComputedContacts : Result container returned by compute()
    """

    def __init__(
        self,
        query,
        database,
        params=None,
        compute_strategy: Literal["default"] = "default",
    ):
        self.query = query
        self.database = database

        self._contact_computers = {"default": SerialContacts}
        self._compute_strategy = compute_strategy
        self._contact_strategy = ExactContacts

        self.params = params if params is not None else DEFAULT_SIM_PARAMS

    def compute(
        self, strategy_or_computer=None, start=None, stop=None, step=1, **kwargs
    ):
        """Compute contacts between query and database.

        Parameters
        ----------
        strategy_or_computer : str, optional
            Computation strategy name. Default: "default".
        start : int, optional
            First frame to process.
        stop : int, optional
            Last frame to process (exclusive).
        step : int, default=1
            Frame step size.
        **kwargs : dict
            Additional arguments passed to contact computer (e.g., cutoff).

        Returns
        -------
        ComputedContacts
            Container with contact data and analysis methods.

        Raises
        ------
        ValueError
            If strategy_or_computer is not recognized.
        """
        if strategy_or_computer is None:
            strategy_or_computer = self._compute_strategy

        contact_computer_class = self._contact_computers.get(strategy_or_computer, None)
        if contact_computer_class is None:
            strats = ", ".join(self._contact_computers.keys())
            raise ValueError(
                f"Unknown strategy or computer: {strategy_or_computer}. Available strategies are: {strats}."
            )

        logger.debug("Using contact computer strategy: %s", strategy_or_computer)

        contact_computer = contact_computer_class(
            self.query.universe, self.query, self.database, **kwargs
        )
        contact_computer.run(verbose=True, start=start, stop=stop, step=step)
        contact_frames = contact_computer.contact_frames

        logger.info("Aggregating contacts into durations...")
        contact_strategy_instance = self._contact_strategy(
            self.query.universe,
            contact_frames,
            self.params.get("norm_factor", 1.0),
        )
        contact_strategy_instance.run()
        logger.info("Contact aggregation complete")

        return ComputedContacts(contact_strategy_instance, self)
