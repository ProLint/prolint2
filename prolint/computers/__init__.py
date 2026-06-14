"""Contact computation algorithms.

This module provides classes for computing distance-based contacts
between atom groups in molecular dynamics trajectories.
"""

from prolint.computers.contacts import SerialContacts
from prolint.computers.base import ContactComputerBase

__all__ = ["SerialContacts", "ContactComputerBase"]
