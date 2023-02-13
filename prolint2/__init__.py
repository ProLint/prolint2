"""ProLint2: Lipid-Protein Interaction Analysis"""

# Add imports here
import os
from .prolint2 import *
from .interactive_sel import *
from .visualization import *

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# to get the paths relative to the root of the package
_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data():
    return os.path.join(_ROOT, "data")


# to get the path to the config file
def get_config():
    return os.path.join(_ROOT, "config.ini")
