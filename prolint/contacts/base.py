"""Base contact storage module.

This module provides the abstract base class for contact storage strategies.
"""

from collections import defaultdict

from typing import List, Union, Callable


class BaseContactStore:
    """Abstract base class for contact storage strategies.

    Defines the interface for storing and computing metrics on contact data.
    Subclasses implement specific storage and aggregation strategies.

    Parameters
    ----------
    ts : Universe
        MDAnalysis Universe instance.
    contact_frames : dict
        Nested dict mapping residue_id -> database_id -> list of frame indices.
    norm_factor : float, default=1.0
        Normalization factor for duration calculations.

    See Also
    --------
    ExactContacts : Concrete implementation for exact contact aggregation
    """

    def __init__(self, ts, contact_frames, norm_factor: float = 1.0):
        self.norm_factor = float(norm_factor)
        self.contact_frames = contact_frames

        self._universe = ts
        self._contacts = defaultdict(lambda: defaultdict(dict))

    def run(self, database_resnames: Union[str, List] = None):
        """Process contact frames into aggregated contact data.

        Parameters
        ----------
        database_resnames : str or list of str, optional
            Residue names to process. If None, processes all.

        Raises
        ------
        NotImplementedError
            Subclasses must implement this method.
        """
        raise NotImplementedError("Subclasses should implement this method.")

    def compute(self, metric: str, target_resname=None):
        """Compute a standard metric on contacts.

        Parameters
        ----------
        metric : {"max", "sum", "mean", "occupancy"}
            Metric to compute.
        target_resname : str, optional
            Filter by residue name.

        Returns
        -------
        dict
            Computed metric values.

        Raises
        ------
        ValueError
            If metric is not recognized.
        """
        if metric in ["max", "sum", "mean", "occupancy"]:
            return self.compute_metric(metric, target_resname)
        else:
            raise ValueError(
                "Invalid metric specified. Use 'max', 'sum', 'mean', 'occupancy. For more complex metrics, use `apply_function`."
            )

    def compute_metric(self, metric: str, target_resname=None):
        """Compute a specific metric. Implemented by subclasses.

        Parameters
        ----------
        metric : str
            Metric to compute.
        target_resname : str, optional
            Filter by residue name.

        Raises
        ------
        NotImplementedError
            Subclasses must implement this method.
        """
        raise NotImplementedError("Subclasses should implement this method.")

    def apply_function(self, func: Callable, target_resname=None):
        """Apply a custom function to contact data.

        Parameters
        ----------
        func : callable
            Function to apply to contact durations.
        target_resname : str, optional
            Filter by residue name.

        Raises
        ------
        NotImplementedError
            Subclasses must implement this method.
        """
        raise NotImplementedError("Subclasses should implement this method.")

    @property
    def contacts(self):
        """Processed contact data.

        Returns
        -------
        dict
            Contact data organized by residue and database molecule.

        Raises
        ------
        ValueError
            If run() has not been called yet.
        """
        if not self._contacts:
            raise ValueError("No contacts have been computed yet. Call run() first.")
        return self._contacts
