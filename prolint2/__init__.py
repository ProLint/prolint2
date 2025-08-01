"""ProLint2: Lipid-Protein Interaction Analysis"""

# Add imports here
import os

# Handle versioneer
from ._version import get_versions
from .core import Universe

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
