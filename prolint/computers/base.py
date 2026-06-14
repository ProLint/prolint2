"""Base contact computation module.

This module provides the abstract base class for contact detection algorithms.
"""

from abc import ABC
from MDAnalysis.analysis.base import AnalysisBase


class ContactComputerBase(AnalysisBase, ABC):
    """Abstract base class for contact computation algorithms.

    Extends MDAnalysis AnalysisBase with contact-specific operations.
    Subclasses implement specific algorithms for detecting contacts
    between atom groups.

    See Also
    --------
    SerialContacts : Concrete implementation using grid-based search
    """

    def __add__(self, other):
        """Combine contacts from two computations."""
        pass

    def intersection(self, other):
        """Find contacts common to both computations."""
        pass

    def union(self, other):
        """Combine all contacts from both computations."""
        pass
