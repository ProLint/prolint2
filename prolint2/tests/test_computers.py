import warnings
import pytest
import numpy as np
from collections import defaultdict

import MDAnalysis as mda
from prolint2 import Universe
from prolint2.sampledata import GIRKDataSample, COX1DataSample, SMODataSample
from prolint2.computers.contacts import SerialContacts

GIRK = GIRKDataSample()
COX1 = COX1DataSample()
SMO = SMODataSample()

warnings.filterwarnings("ignore")


@pytest.fixture(scope="module")
def universes():
    return [
        Universe(GIRK.coordinates, GIRK.trajectory),
        Universe(COX1.coordinates, COX1.trajectory),
        Universe(SMO.coordinates, SMO.trajectory),
    ]


class TestSerialContacts:
    def test_init(self, universes):
        for u in universes:
            cutoff = 6.0
            contacts = SerialContacts(u, u.query, u.database, cutoff)
            assert contacts.query == u.query
            assert contacts.database == u.database
            assert contacts.cutoff == cutoff
            assert len(contacts.q_resids) == len(u.query.resids)
            assert len(contacts.db_resids) == len(u.database.resids)
            assert len(contacts.db_resnames) == len(u.database.resnames)
            assert contacts.contacts is None
            assert isinstance(contacts.contact_frames, defaultdict)

    def test_init_invalid_selection(self):
        with pytest.raises(mda.exceptions.SelectionError):
            universe = Universe(GIRK.coordinates, GIRK.trajectory)
            query = universe.select_atoms("not protein")
            database = universe.select_atoms("lipid")
            SerialContacts(universe, query, database)

    def test_init_invalid_cutoff(self):
        with pytest.raises(ValueError):
            universe = Universe(GIRK.coordinates, GIRK.trajectory)
            cutoff = -2.0
            SerialContacts(universe, universe.query, universe.database, cutoff)

    def test_validate_inputs(self, universes):
        for u in universes:
            contacts = SerialContacts(u, u.query, u.database)
            contacts._validate_inputs()  # No exception should be raised

            query_empty = u.select_atoms("protein and name XYZ")
            with pytest.raises(ValueError):
                SerialContacts(u, query_empty, u.database)

            cutoff_zero = 0.0
            with pytest.raises(ValueError):
                SerialContacts(u, u.query, u.database, cutoff_zero)

    def test_get_residue_lipid_info(self, universes):
        for u in universes:
            contacts = SerialContacts(u, u.query, u.database)
            pair = (0, 1)
            residue_id, lipid_id, lipid_name = contacts._get_residue_lipid_info(pair)
            assert isinstance(residue_id, np.int64)
            assert isinstance(lipid_id, np.int64)
            assert isinstance(lipid_name, str)

    def test_compute_pairs(self, universes):
        for u in universes:
            contacts = SerialContacts(u, u.query, u.database)
            pairs = contacts._compute_pairs()
            assert isinstance(pairs, np.ndarray)
            assert pairs.shape[1] == 2
