"""Unit tests for prolint.core.groups module."""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from abc import ABC


class TestPLAtomGroupBase:
    """Tests for PLAtomGroupBase abstract class."""

    def test_is_abstract_base_class(self):
        """PLAtomGroupBase is an ABC with required abstract methods."""
        from prolint.core.groups import PLAtomGroupBase
        import inspect

        assert issubclass(PLAtomGroupBase, ABC)

        abstract_methods = {
            name for name, method in inspect.getmembers(PLAtomGroupBase)
            if getattr(method, "__isabstractmethod__", False)
        }
        expected = {"add", "remove", "get_resnames", "get_resids", "get_all_resids",
                    "filter_resids_by_resname", "unique_resnames", "resname_counts"}
        assert abstract_methods == expected


class TestExtendedAtomGroup:
    """Tests for ExtendedAtomGroup class."""

    def test_inherits_from_mda_atomgroup_and_base(self):
        """ExtendedAtomGroup inherits from MDAnalysis AtomGroup and PLAtomGroupBase."""
        from prolint.core.groups import ExtendedAtomGroup, PLAtomGroupBase
        import MDAnalysis as mda

        assert issubclass(ExtendedAtomGroup, mda.AtomGroup)
        assert issubclass(ExtendedAtomGroup, PLAtomGroupBase)

    def test_build_selection_string(self):
        """_build_selection_string creates correct MDAnalysis selection strings."""
        from prolint.core.groups import ExtendedAtomGroup

        mock = MagicMock(spec=ExtendedAtomGroup)
        mock._build_selection_string = lambda **kw: ExtendedAtomGroup._build_selection_string(mock, **kw)

        assert "resname CHOL" in mock._build_selection_string(resname="CHOL")
        assert "name CA" in mock._build_selection_string(atomname="CA")
        assert "resid 42" in mock._build_selection_string(resnum=42)

        with pytest.raises(ValueError):
            mock._build_selection_string()  # No criteria

    def test_get_resids_filters_by_resname(self):
        """get_resids returns residue IDs matching a resname."""
        from prolint.core.groups import ExtendedAtomGroup

        mock = MagicMock(spec=ExtendedAtomGroup)
        mock.residues.resids = np.array([1, 2, 3, 4, 5])
        mock.residues.resnames = np.array(["ALA", "GLY", "ALA", "VAL", "ALA"])

        result = ExtendedAtomGroup.get_resids(mock, "ALA")
        np.testing.assert_array_equal(result, np.array([1, 3, 5]))

    def test_filter_resids_by_resname(self):
        """filter_resids_by_resname filters residue IDs by resname."""
        from prolint.core.groups import ExtendedAtomGroup

        mock = MagicMock(spec=ExtendedAtomGroup)
        mock._stored_resnames = np.array(["ALA", "GLY", "ALA"])
        mock._stored_resids = np.array([1, 2, 3])

        result = ExtendedAtomGroup.filter_resids_by_resname(mock, [1, 2, 3], "ALA")
        np.testing.assert_array_equal(result, np.array([1, 3]))
