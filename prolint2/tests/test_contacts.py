import numpy as np
import pytest
import MDAnalysis as mda
from prolint2.sampledata import GIRK
import re
from prolint2 import QueryProteins, MembraneDatabase
from prolint2.contacts import SerialContacts, SerialDistances, Contacts

@pytest.fixture
def universe():
    return mda.Universe(GIRK.coordinates, GIRK.trajectory)

def test_init(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    cutoff = 8.0

    contacts = SerialContacts(universe, query, database, cutoff)
    assert contacts.query == query
    assert contacts.database == database
    assert contacts.cutoff == cutoff
    assert contacts.q_resids == query.resindices.tolist()
    assert contacts.db_resids == database.resindices.tolist()
    assert (contacts.db_resnames == database.resnames).all()
    assert np.array_equal(contacts.dp_resnames_unique, np.unique(database.resnames))

def test_init_exception(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("")
    cutoff = 8.0
    
    with pytest.raises(ValueError, match=re.escape("Invalid selection. Empty AtomGroup(s).")):
        contacts = SerialContacts(universe, query, database, cutoff)
    
    query = universe.select_atoms("")
    database = universe.select_atoms("resname CHOL")
    cutoff = 8.0
    
    with pytest.raises(ValueError, match=re.escape("Invalid selection. Empty AtomGroup(s).")):
        contacts = SerialContacts(universe, query, database, cutoff)
        
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    cutoff = -8.0
    
    with pytest.raises(ValueError, match=re.escape("The cutoff must be greater than 0.")):
        contacts = SerialContacts(universe, query, database, cutoff)
        
def test_prepare(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    cutoff = 8.0

    contacts = SerialContacts(universe, query, database, cutoff)
    contacts._prepare()
    q_resids = contacts.q_resids
    dp_resnames_unique = contacts.dp_resnames_unique
    assert contacts.contacts == {k: {v: [] for v in dp_resnames_unique} for k in q_resids}
    assert contacts.contacts_sum == {k: {v: 0 for v in dp_resnames_unique} for k in q_resids}
    assert contacts.contact_frames == {}

def test_single_frame(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    cutoff = 8.0

    sc = SerialContacts(universe, query, database, cutoff)
    sc._prepare()
    # sc._single_frame()

    # Test if the right number of contacts is computed
    # assert len(sc.contacts) == len(query)
    for k in sc.contacts:
        assert len(sc.contacts[k]) == len(np.unique(database.resnames))

    # Test if the contact_frames attribute is correctly computed
    assert len(sc.contact_frames) <= len(query) * len(database)
    for k, v in sc.contact_frames.items():
        residue_id, lipid_id = map(int, k.split(','))
        lipid_name = database[lipid_id].resname
        assert residue_id in sc.contacts
        assert lipid_name in sc.contacts[residue_id]
        assert lipid_id in sc.contacts[residue_id][lipid_name]
        assert len(v) <= len(universe.universe.trajectory)
        for frame in v:
            assert frame in universe.universe.trajectory.frame_indices

def test_SerialDistances_instance(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    lipid_id = 0
    residue_id = 0
    frame_filter = np.array([0, 1, 2])
    sd = SerialDistances(universe, query, database, lipid_id, residue_id, frame_filter)
    assert isinstance(sd, mda.analysis.base.AnalysisBase)

def test_SerialDistances_empty_atomgroup_selection(universe):
    query = universe.select_atoms("protein")
    database = universe.select_atoms("")
    lipid_id = 0
    residue_id = 0
    frame_filter = np.array([0, 1, 2])
    with pytest.raises(ValueError, match="Invalid selection. Empty AtomGroup\(s\)."):
        sd = SerialDistances(universe, query, database, lipid_id, residue_id, frame_filter)

def test_SerialDistances_prepare(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    lipid_id = 0
    residue_id = 0
    frame_filter = np.array([0, 1, 2])
    sd = SerialDistances(universe, query, database, lipid_id, residue_id, frame_filter)
    sd._prepare()
    assert isinstance(sd.result_array, np.ndarray)
    assert sd.result_array.shape == (3, len(sd.lipid_atomgroup), len(sd.resid_atomgroup))

def test_SerialDistances_single_frame(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    lipid_id = 0
    residue_id = 0
    frame_filter = np.array([0, 1, 2])
    sd = SerialDistances(universe, query, database, lipid_id, residue_id, frame_filter)
    sd._prepare()
    # sd._single_frame()
    assert isinstance(sd.result_array[0], np.ndarray)

def test_init_method_sets_correct_attributes(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    query = QueryProteins(query)
    database = MembraneDatabase(database)
    contacts = Contacts(query, database)

    assert contacts.query == query
    assert contacts.database == database
    assert contacts.cutoff is None
    assert contacts.contacts is None

def test_compute_method_sets_cutoff_attribute(universe):
    query = universe.select_atoms("resname LYS")
    database = universe.select_atoms("resname CHOL")
    query = QueryProteins(query)
    database = MembraneDatabase(database)
    contacts = Contacts(query, database)
    contacts.compute(7)

    assert contacts.cutoff == 7

