import warnings
import pandas as pd
import pytest

import MDAnalysis as mda
from typing import Literal, get_args

from prolint2 import Universe
from prolint2.core.groups import ExtendedAtomGroup
from prolint2.sampledata import GIRKDataSample, COX1DataSample, SMODataSample

GIRK = GIRKDataSample()
COX1 = COX1DataSample()
SMO = SMODataSample()

warnings.filterwarnings("ignore")

TimeUnitLiteral = Literal["fs", "ps", "ns", "us", "ms", "s"]

# Build VALID_UNITS from TimeUnitLiteral
VALID_UNITS = get_args(TimeUnitLiteral)


@pytest.fixture(scope="module")
def universes():
    return [
        Universe(GIRK.coordinates, GIRK.trajectory),
        Universe(COX1.coordinates, COX1.trajectory),
        Universe(SMO.coordinates, SMO.trajectory),
    ]


class TestUniverse:
    def test_init_with_universe(self, universes):
        for u in universes:
            assert isinstance(u, Universe)
            assert isinstance(u, mda.Universe)

    def test_query_property(self, universes):
        for u in universes:
            assert u.query is not None
            assert isinstance(u.query, ExtendedAtomGroup)

    def test_database_property(self, universes):
        for u in universes:
            assert u.database is not None
            assert isinstance(u.database, ExtendedAtomGroup)

    # def test_compute_contacts(self, universes):
    #     dict_of_contacts = {0: GIRK, 1: COX1, 2: SMO}
    #     for i, u in enumerate(universes):
    #         contacts = u.compute_contacts(cutoff=7)
    #         calc_df = contacts.create_dataframe(u.trajectory.n_frames)
    #         ref_df = pd.read_csv(dict_of_contacts[i].contacts, index_col=[0, 1])
    #         ref_df.columns = ref_df.columns.astype(int)
    #         ref_df = ref_df.astype("int8")
    #         pd.testing.assert_frame_equal(calc_df, ref_df)

    # def test_load_contacts_from_file(self, universes):
    #     dict_of_contacts = {0: GIRK, 1: COX1, 2: SMO}
    #     for i, u in enumerate(universes):
    #         contacts = u.load_contacts_from_file(dict_of_contacts[i].contacts)
    #         calc_df = contacts.create_dataframe(u.trajectory.n_frames)
    #         ref_df = pd.read_csv(dict_of_contacts[i].contacts, index_col=[0, 1])
    #         ref_df.columns = ref_df.columns.astype(int)
    #         ref_df = ref_df.astype("int8")
    #         pd.testing.assert_frame_equal(calc_df, ref_df)

    # def test_units_property(self, universes):
    #     for u in universes:
    #         assert u.units in VALID_UNITS

    # def test_units_property_setter(self, universes):
    #     for u in universes:
    #         new_units = "ns"
    #         u.units = new_units
    #         assert u.units == new_units
    #         assert u.units in VALID_UNITS

    # def test_normalize_by_property(self, universes):
    #     for u in universes:
    #         assert u.normalize_by in ["counts", "actual_time", "time_fraction"]

    # def test_normalize_by_property_setter(self, universes):
    #     for u in universes:
    #         new_normalize_by = "counts"
    #         u.normalize_by = new_normalize_by
    #         assert u.normalize_by == new_normalize_by
    #         assert u.normalize_by in ["counts", "actual_time", "time_fraction"]

    # def test_query_property_setter(self, universes):
    #     for u in universes:
    #         new_query = u.select_atoms("protein and (name CA or name BB)")
    #         u.query = new_query
    #         assert u.query is not None
    #         assert isinstance(u.query, ExtendedAtomGroup)

    # def test_database_property_setter(self, universes):
    #     for u in universes:
    #         new_database = u.select_atoms("resname CHOL or resname CHL1")
    #         u.database = new_database
    #         assert u.database is not None
    #         assert isinstance(u.database, ExtendedAtomGroup)
