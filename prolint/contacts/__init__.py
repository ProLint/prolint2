"""Contact storage and aggregation strategies.

This module provides classes for storing contact data and computing
duration-based metrics from contact frame indices.
"""

from prolint.contacts.base import BaseContactStore
from prolint.contacts.exact_contacts import ExactContacts

__all__ = ["BaseContactStore", "ExactContacts"]
