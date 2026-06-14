"""Unit tests for prolint.contacts.exact_contacts module."""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch


class TestExactContacts:
    """Tests for ExactContacts class."""

    def test_inherits_from_base_contact_store(self):
        """ExactContacts inherits from BaseContactStore."""
        from prolint.contacts.exact_contacts import ExactContacts
        from prolint.contacts.base import BaseContactStore

        assert issubclass(ExactContacts, BaseContactStore)

    def test_run_processes_database_resnames(self):
        """run() processes database residue names."""
        from prolint.contacts.exact_contacts import ExactContacts

        mock_ts = MagicMock()
        mock_ts.database.filter_resids_by_resname.return_value = np.array([10])
        contact_frames = {1: {10: [0, 1, 2]}}

        store = ExactContacts(mock_ts, contact_frames)

        with patch.object(store, "compute_database_durations", return_value={10: np.array([3])}):
            store.run(database_resnames=["CHOL"])

        assert 1 in store._contacts

    def test_compute_metric_returns_correct_values(self):
        """compute_metric returns correct metric values."""
        from prolint.contacts.exact_contacts import ExactContacts

        mock_ts = MagicMock()
        mock_ts.trajectory.n_frames = 100

        store = ExactContacts(mock_ts, {})
        store._contacts[1]["CHOL"] = {10: np.array([3, 3]), 11: np.array([2, 2])}

        result = store.compute_metric("sum", target_resname="CHOL")

        assert result[1]["CHOL"]["per_id"][10] == 6  # 3 + 3
        assert result[1]["CHOL"]["per_id"][11] == 4  # 2 + 2

    def test_compute_metric_raises_for_unknown(self):
        """compute_metric raises ValueError for unknown metric."""
        from prolint.contacts.exact_contacts import ExactContacts

        store = ExactContacts(MagicMock(), {})
        store._contacts[1]["CHOL"] = {10: np.array([3])}

        with pytest.raises(ValueError, match="Unknown metric"):
            store.compute_metric("invalid")

    def test_apply_function_applies_custom_function(self):
        """apply_function applies custom function to durations."""
        from prolint.contacts.exact_contacts import ExactContacts

        store = ExactContacts(MagicMock(), {})
        store._contacts[1]["CHOL"] = {10: np.array([3, 5, 2])}

        result = store.apply_function(len)

        assert result[1]["CHOL"][10] == 3  # 3 binding events
