"""Unit tests for prolint.computers.contacts module."""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from collections import defaultdict


class TestSerialContacts:
    """Tests for SerialContacts class."""

    def test_inherits_from_contact_computer_base(self):
        """SerialContacts inherits from ContactComputerBase."""
        from prolint.computers.contacts import SerialContacts
        from prolint.computers.base import ContactComputerBase

        assert issubclass(SerialContacts, ContactComputerBase)

    def test_validates_empty_query_and_database(self):
        """_validate_inputs raises for empty AtomGroups."""
        from prolint.computers.contacts import SerialContacts

        sc = SerialContacts.__new__(SerialContacts)
        sc.query = MagicMock()
        sc.query.__len__ = MagicMock(return_value=0)
        sc.database = MagicMock()
        sc.database.__len__ = MagicMock(return_value=100)
        sc.cutoff = 7.0

        with pytest.raises(ValueError, match="Empty AtomGroup"):
            sc._validate_inputs()

    def test_validates_cutoff(self):
        """_validate_inputs raises for invalid cutoff."""
        from prolint.computers.contacts import SerialContacts

        sc = SerialContacts.__new__(SerialContacts)
        sc.query = MagicMock()
        sc.query.__len__ = MagicMock(return_value=100)
        sc.database = MagicMock()
        sc.database.__len__ = MagicMock(return_value=100)
        sc.cutoff = 0

        with pytest.raises(ValueError, match="cutoff must be greater than 0"):
            sc._validate_inputs()

    def test_single_frame_stores_contacts(self):
        """_single_frame stores contact pairs in contact_frames."""
        from prolint.computers.contacts import SerialContacts

        sc = SerialContacts.__new__(SerialContacts)
        sc.contact_frames = defaultdict(lambda: defaultdict(list))
        sc._last_log_frame = 0
        sc._log_interval = 10
        sc.n_frames = 100
        sc._frame_index = 0
        sc._ts = MagicMock()
        sc._ts.frame = 0

        sc.query = MagicMock()
        sc.query.resids = np.array([1, 2, 3])
        sc.database = MagicMock()
        sc.database.resids = np.array([10, 11, 12])

        pairs = np.array([[0, 0], [1, 2]])
        with patch.object(sc, "_compute_pairs", return_value=pairs):
            with patch("prolint.computers.contacts.fast_unique_comparison") as mock_unique:
                mock_unique.return_value = (
                    np.array([1, 2]),
                    np.array([10, 12]),
                    np.array(["POPC", "CHOL"]),
                )
                sc._single_frame()

        assert 0 in sc.contact_frames[1][10]
        assert 0 in sc.contact_frames[2][12]
