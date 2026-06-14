"""Analysis base module.

This module provides the base classes and registry for ProLint analyses.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Type, List, TYPE_CHECKING

if TYPE_CHECKING:
    from prolint.core.contact_provider import ComputedContacts


@dataclass
class AnalysisResult:
    """Container for analysis results.

    The ``data`` attribute holds primary analysis data (varies by analysis type).
    The ``metadata`` attribute holds metadata about the analysis (parameters, timestamps, etc.).
    """

    data: Dict[str, Any] = field(default_factory=dict)
    """Primary analysis data (varies by analysis type)."""

    metadata: Dict[str, Any] = field(default_factory=dict)
    """Metadata about the analysis (parameters, timestamps, etc.)."""


class BaseAnalysis(ABC):
    """Abstract base class for all ProLint analyses.

    Provides common functionality for filtering contacts, building
    residue mappings, and generating frame ranges.

    Parameters
    ----------
    universe : Universe
        ProLint Universe instance.
    contacts : ComputedContacts
        Computed contact data to analyze.

    See Also
    --------
    AnalysisRegistry : Registry for creating analyses by name
    ComputedContacts.analyze : High-level interface to run analyses
    """

    name: str = "base_analysis"
    """Analysis name for registry."""

    description: str = "Base analysis class"
    """Human-readable description."""

    def __init__(self, universe, contacts: "ComputedContacts"):
        self.universe = universe
        self.contacts = contacts
        self._validate_inputs()

    def _validate_inputs(self):
        if self.universe is None:
            raise ValueError("Universe cannot be None")
        if self.contacts is None:
            raise ValueError("contacts cannot be None")

    @abstractmethod
    def run(self, **kwargs) -> AnalysisResult:
        """Run the analysis and return results.

        Parameters
        ----------
        **kwargs : dict
            Analysis-specific parameters.

        Returns
        -------
        AnalysisResult
            Container with analysis data and metadata.
        """
        pass

    def _get_database_id_to_resname(self) -> Dict[int, str]:
        """Build mapping from database residue ID to residue name.

        Returns
        -------
        dict
            Mapping of database residue ID to residue name string.
        """
        return self.universe.database.get_resnames(
            self.universe.database.residues.resids, out=dict
        )

    def _filter_by_database_type(
        self, database_type: Optional[str] = None
    ) -> Dict[int, Dict[int, List[int]]]:
        """Filter contact_frames to only include specified database type.

        Parameters
        ----------
        database_type : str, optional
            Database residue name to filter by (e.g., "CHOL").
            If None, returns all contact frames unfiltered.

        Returns
        -------
        dict
            Nested dict mapping query_resid -> database_id -> list of frame indices.
        """
        if database_type is None:
            return self.contacts.contact_frames

        id_to_resname = self._get_database_id_to_resname()
        filtered = {}

        for query_resid, db_dict in self.contacts.contact_frames.items():
            filtered_db = {}
            for db_id, frames in db_dict.items():
                if id_to_resname.get(db_id) == database_type:
                    filtered_db[db_id] = frames
            if filtered_db:
                filtered[query_resid] = filtered_db

        return filtered

    def _get_frame_range(
        self,
        frame_start: Optional[int] = None,
        frame_end: Optional[int] = None,
        frame_step: int = 1,
    ) -> List[int]:
        """Generate list of frame indices for analysis.

        Parameters
        ----------
        frame_start : int, optional
            First frame index. Defaults to 0.
        frame_end : int, optional
            Last frame index (exclusive). Defaults to total frames.
        frame_step : int, default=1
            Step between frames.

        Returns
        -------
        list of int
            Frame indices to process.
        """
        n_frames = self.universe.trajectory.n_frames
        start = frame_start if frame_start is not None else 0
        end = frame_end if frame_end is not None else n_frames
        return list(range(start, end, frame_step))


class AnalysisRegistry:
    """Registry for analysis types.

    Manages registration and creation of analysis classes. All built-in
    analyses are registered automatically on import.

    Examples
    --------
    List available analyses:

    >>> from prolint.analysis import AnalysisRegistry
    >>> AnalysisRegistry.available()
    ['timeseries', 'database_contacts', 'kinetics', ...]

    Create an analysis:

    >>> analysis = AnalysisRegistry.create("timeseries", universe, contacts)
    >>> result = analysis.run(database_type="CHOL")
    """

    _registry: Dict[str, Type[BaseAnalysis]] = {}

    @classmethod
    def register(cls, name: str, analysis_class: Type[BaseAnalysis]):
        """Register an analysis class.

        Parameters
        ----------
        name : str
            Name to register under.
        analysis_class : type
            Analysis class (must inherit from BaseAnalysis).
        """
        if not issubclass(analysis_class, BaseAnalysis):
            raise TypeError(f"{analysis_class} must be a subclass of BaseAnalysis")
        cls._registry[name] = analysis_class

    @classmethod
    def create(
        cls, name: str, universe, contacts: "ComputedContacts", **kwargs
    ) -> BaseAnalysis:
        """Create an analysis instance.

        Parameters
        ----------
        name : str
            Analysis type name.
        universe : Universe
            ProLint Universe instance.
        contacts : ComputedContacts
            Computed contact data.
        **kwargs : dict
            Additional arguments for the analysis.

        Returns
        -------
        BaseAnalysis
            Initialized analysis instance.
        """
        if name not in cls._registry:
            available = ", ".join(cls._registry.keys())
            raise ValueError(f"Unknown analysis: {name}. Available: {available}")
        return cls._registry[name](universe, contacts, **kwargs)

    @classmethod
    def available(cls) -> List[str]:
        """List available analysis types.

        Returns
        -------
        list of str
            Registered analysis names.
        """
        return list(cls._registry.keys())

    @classmethod
    def get_class(cls, name: str) -> Type[BaseAnalysis]:
        """Get an analysis class by name.

        Parameters
        ----------
        name : str
            Analysis name.

        Returns
        -------
        type
            Analysis class.
        """
        if name not in cls._registry:
            raise ValueError(f"Unknown analysis: {name}")
        return cls._registry[name]
