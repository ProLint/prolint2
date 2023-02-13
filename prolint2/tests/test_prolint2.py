"""
Unit tests for the prolint2 package.
"""

# Import package, test suite, and other packages as needed
import sys

from pandas.util.testing import assert_frame_equal

import prolint2
from prolint2.sampledata import GIRK


def test_prolint2_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "prolint2" in sys.modules


def test_prolint2_initialize():
    """Initialize test, will always pass so long as the prolint2.PL2 object can be initialized."""
    assert prolint2.PL2(GIRK.coordinates, GIRK.trajectory)
