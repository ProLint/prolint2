"""Unit tests for prolint.core.contact_provider module."""

import pytest
from unittest.mock import MagicMock, patch


class TestComputedContacts:
    """Tests for ComputedContacts class."""

    def test_properties_return_stored_values(self):
        """ComputedContacts properties return stored values."""
        from prolint.core.contact_provider import ComputedContacts

        mock_strategy = MagicMock()
        mock_strategy.contacts = {"resid1": {"CHOL": {"mol1": [1, 2, 3]}}}
        mock_strategy.contact_frames = {1: {10: [0, 1, 2]}}
        mock_strategy.norm_factor = 2.5
        mock_provider = MagicMock()

        cc = ComputedContacts(mock_strategy, mock_provider)

        assert cc.contacts == mock_strategy.contacts
        assert cc.contact_frames == {1: {10: [0, 1, 2]}}
        assert cc.norm_factor == 2.5

    def test_compute_metric_delegates_to_strategy(self):
        """compute_metric delegates to contact strategy."""
        from prolint.core.contact_provider import ComputedContacts

        mock_strategy = MagicMock()
        mock_strategy.compute.return_value = {1: {"CHOL": {"global": 0.5}}}
        mock_provider = MagicMock()

        cc = ComputedContacts(mock_strategy, mock_provider)
        result = cc.compute_metric("max", target_resname="CHOL")

        mock_strategy.compute.assert_called_once_with("max", target_resname="CHOL")

    def test_analyze_uses_registry(self):
        """analyze uses AnalysisRegistry to create analysis."""
        from prolint.core.contact_provider import ComputedContacts

        mock_strategy = MagicMock()
        mock_provider = MagicMock()
        mock_provider.query.universe = MagicMock()

        cc = ComputedContacts(mock_strategy, mock_provider)

        with patch("prolint.analysis.AnalysisRegistry") as MockRegistry:
            mock_analysis = MagicMock()
            MockRegistry.create.return_value = mock_analysis

            cc.analyze("timeseries", database_type="CHOL")

            MockRegistry.create.assert_called_once()
            mock_analysis.run.assert_called_once()


class TestContactsProvider:
    """Tests for ContactsProvider class."""

    def test_stores_query_and_database(self):
        """ContactsProvider stores query and database references."""
        from prolint.core.contact_provider import ContactsProvider

        mock_query = MagicMock()
        mock_database = MagicMock()

        provider = ContactsProvider(mock_query, mock_database)

        assert provider.query is mock_query
        assert provider.database is mock_database

    def test_compute_returns_computed_contacts(self):
        """compute returns ComputedContacts instance."""
        from prolint.core.contact_provider import ContactsProvider, ComputedContacts

        mock_query = MagicMock()
        mock_query.universe = MagicMock()
        mock_database = MagicMock()

        provider = ContactsProvider(mock_query, mock_database)

        with patch.object(provider, "_contact_computers") as mock_computers:
            mock_computer_cls = MagicMock()
            mock_computer = MagicMock()
            mock_computer.contact_frames = {1: {10: [0, 1, 2]}}
            mock_computer_cls.return_value = mock_computer
            mock_computers.get.return_value = mock_computer_cls

            with patch.object(provider, "_contact_strategy") as mock_strategy:
                mock_strategy.return_value = MagicMock()
                result = provider.compute(cutoff=7.0)

                assert isinstance(result, ComputedContacts)

    def test_compute_raises_for_unknown_strategy(self):
        """compute raises ValueError for unknown strategy."""
        from prolint.core.contact_provider import ContactsProvider

        provider = ContactsProvider(MagicMock(), MagicMock())

        with pytest.raises(ValueError, match="Unknown strategy"):
            provider.compute(strategy_or_computer="nonexistent")
