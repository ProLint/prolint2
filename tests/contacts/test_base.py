"""Unit tests for prolint.contacts.base module."""

import pytest
from unittest.mock import MagicMock, patch


class TestBaseContactStore:
    """Tests for BaseContactStore class."""

    def test_stores_initialization_params(self):
        """BaseContactStore stores universe, contact_frames, and norm_factor."""
        from prolint.contacts.base import BaseContactStore

        mock_ts = MagicMock()
        contact_frames = {1: {10: [0, 1, 2]}}

        store = BaseContactStore(mock_ts, contact_frames, norm_factor=2.5)

        assert store._universe is mock_ts
        assert store.contact_frames is contact_frames
        assert store.norm_factor == 2.5

    def test_run_raises_not_implemented(self):
        """run() raises NotImplementedError."""
        from prolint.contacts.base import BaseContactStore

        store = BaseContactStore(MagicMock(), {})

        with pytest.raises(NotImplementedError):
            store.run()

    def test_compute_delegates_valid_metrics(self):
        """compute() delegates to compute_metric for valid metrics."""
        from prolint.contacts.base import BaseContactStore

        store = BaseContactStore(MagicMock(), {})

        with patch.object(store, "compute_metric", return_value={}) as mock:
            store.compute("max", target_resname="CHOL")
            mock.assert_called_once_with("max", "CHOL")

    def test_compute_rejects_invalid_metrics(self):
        """compute() raises ValueError for invalid metrics."""
        from prolint.contacts.base import BaseContactStore

        store = BaseContactStore(MagicMock(), {})

        with pytest.raises(ValueError, match="Invalid metric"):
            store.compute("invalid_metric")

    def test_contacts_property_raises_when_empty(self):
        """contacts property raises ValueError when no contacts computed."""
        from prolint.contacts.base import BaseContactStore

        store = BaseContactStore(MagicMock(), {})

        with pytest.raises(ValueError, match="No contacts"):
            _ = store.contacts
