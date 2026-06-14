"""Core module for ProLint.

Provides the main classes for biomolecular interaction analysis:
- Universe: Entry point extending MDAnalysis Universe
- ExtendedAtomGroup: Enhanced atom group with additional properties
- ContactsProvider: Orchestrates contact computation
"""

from .universe import Universe
from .groups import ExtendedAtomGroup
from .contact_provider import ContactsProvider
