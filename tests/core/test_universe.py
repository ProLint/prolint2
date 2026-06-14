"""Unit tests for prolint.core.universe module."""

import pytest
from unittest.mock import MagicMock, patch


class TestUniverse:
    """Tests for Universe class."""

    def test_inherits_from_mda_universe(self):
        """Universe inherits from MDAnalysis Universe."""
        from prolint.core.universe import Universe
        import MDAnalysis as mda

        assert issubclass(Universe, mda.Universe)

    def test_query_and_database_properties(self, mock_prolint_universe):
        """query and database properties return stored values."""
        assert mock_prolint_universe.query is not None
        assert mock_prolint_universe.database is not None

    def test_setters_reject_invalid_types(self):
        """query/database setters reject non-AtomGroup types."""
        from prolint.core.universe import Universe

        u = Universe.__new__(Universe)
        u._query = MagicMock()
        u._database = MagicMock()

        with pytest.raises(TypeError):
            u.query = "not_an_atom_group"

        with pytest.raises(TypeError):
            u.database = "not_an_atom_group"

    def test_compute_contacts_creates_provider(self):
        """compute_contacts creates ContactsProvider and calls compute."""
        from prolint.core.universe import Universe

        u = Universe.__new__(Universe)
        u._query = MagicMock()
        u._database = MagicMock()
        u.params = {"units": "us", "normalizer": "counts", "norm_factor": 1.0}

        with patch("prolint.core.universe.ContactsProvider") as MockProvider:
            mock_instance = MagicMock()
            MockProvider.return_value = mock_instance

            u.compute_contacts(cutoff=7.0)

            MockProvider.assert_called_once()
            mock_instance.compute.assert_called_once()


class TestValidUnits:
    """Tests for valid units constant."""

    def test_contains_expected_units(self):
        """VALID_UNITS contains expected time units."""
        from prolint.core.universe import VALID_UNITS

        expected = {"fs", "ps", "ns", "us", "ms", "s"}
        assert expected.issubset(set(VALID_UNITS))
