"""
Unit tests for the prolint2 package.
"""

# Import package, test suite, and other packages as needed
import sys
import configparser
import numpy as np
from prolint2 import get_config
from pandas.util.testing import assert_frame_equal
import MDAnalysis as mda
from prolint2 import PL2
from prolint2 import MembraneDatabase
from prolint2 import QueryProteins
from prolint2.sampledata import GIRK
# Getting the config file
config = configparser.ConfigParser(allow_no_value=True)
config.read(get_config())
parameters_config = config["Parameters"]


def test_prolint2_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "prolint2" in sys.modules

def test_prolint2_initialize():
    """Initialize test, will always pass so long as the prolint2.PL2 object can be initialized."""
    assert PL2(GIRK.coordinates, GIRK.trajectory)

def test_init():
    structure = GIRK.coordinates
    trajectory = GIRK.trajectory
    add_lipid_types = []
    pl2 = PL2(structure, trajectory, add_lipid_types)
    
    # Test that the AtomGroup and ResidueGroup have been created
    assert isinstance(pl2.atoms, mda.core.groups.AtomGroup)
    assert isinstance(pl2.residues, mda.core.groups.ResidueGroup)

def test_membrane_residues():
    structure = GIRK.coordinates
    trajectory = GIRK.trajectory
    add_lipid_types = []
    pl2 = PL2(structure, trajectory, add_lipid_types)
    
    # Test that the correct lipid residues have been added to the macros attribute
    lipid_types = parameters_config["lipid_types"].split(", ")
    lipid_types = lipid_types + add_lipid_types
    not_protein_restypes = np.unique(pl2.atoms.select_atoms("not protein").residues.resnames)
    membrane_restypes = []
    for type in lipid_types:
        if type in not_protein_restypes:
            membrane_restypes.append("resname " + type)
    if len(membrane_restypes) == 1:
        membrane_sel = membrane_restypes[0]
    elif len(membrane_restypes) > 1:
        membrane_sel = membrane_restypes[0]
        for type in membrane_restypes[1:]:
            membrane_sel = membrane_sel + " or " + type
    else:
        print("There are not lipid residues in your system")
        
    membrane_residues = pl2.atoms.select_atoms(membrane_sel).residues
    for residue in membrane_residues:
        assert residue.macro == "membrane"

def test_protein_residues():
    structure = GIRK.coordinates
    trajectory = GIRK.trajectory
    add_lipid_types = []
    pl2 = PL2(structure, trajectory, add_lipid_types)
    
    # Test that the correct protein residues have been added to the macros attribute
    protein_sel = "protein"
    if (
        len(pl2.atoms.select_atoms(protein_sel).segments) > 1
        and pl2.atoms.select_atoms(protein_sel).segments.n_atoms   
        == pl2.atoms.select_atoms(protein_sel).n_atoms
    ):
        for segment_idx in range(len(pl2.atoms.select_atoms(protein_sel).segments)):
            protein_segment = pl2.atoms.select_

def test_lipid_types():
    whole = mda.Universe.empty(n_atoms=0, n_residues=0, atom_resindex=([]))
    whole.add_TopologyAttr(
        "resnames", np.array([])
    )
    whole.add_TopologyAttr("resids", np.array([]))
    database = MembraneDatabase(whole)
    lipids = database.lipid_types()
    assert len(lipids) == 0

    # create an AtomGroup with two lipids
    whole = mda.Universe.empty(n_atoms=6, n_residues=2, atom_resindex=([0, 0, 0, 1, 1, 1]))
    whole.add_TopologyAttr(
        "resnames", np.array(["POPC", "POPG"])
    )
    whole.add_TopologyAttr("resids", np.array([0, 1]))
    database = MembraneDatabase(whole)
    lipids = database.lipid_types()
    assert len(lipids) == 2
    assert "POPC" in lipids
    assert "POPG" in lipids

def test_lipid_count():
    # create an AtomGroup with two lipids
    whole = mda.Universe.empty(n_atoms=6, n_residues=2, atom_resindex=([0, 0, 0, 1, 1, 1]))
    whole.add_TopologyAttr(
        "resnames", np.array(["POPC", "POPG"])
    )
    whole.add_TopologyAttr("resids", np.array([0, 1]))
    database = MembraneDatabase(whole)
    lipid_count = database.lipid_count()
    assert len(lipid_count) == 2

def test_str_repr_MembraneDatabase():
    whole = mda.Universe.empty(n_atoms=0)
    database = MembraneDatabase(whole)
    assert str(database) == "<prolint2.MembraneDatabase containing 0 atoms>"
    assert repr(database) == "<prolint2.MembraneDatabase containing 0 atoms>"

    whole = mda.Universe.empty(n_atoms=6)
    database = MembraneDatabase(whole)
    assert str(database) == "<prolint2.MembraneDatabase containing 6 atoms>"
    assert repr(database) == "<prolint2.MembraneDatabase containing 6 atoms>"

def test_str_and_repr_QueryProteins():
    whole = mda.Universe.empty(n_atoms=0)
    database = QueryProteins(whole)
    assert str(database) == "<prolint2.QueryProteins containing 0 atoms>"
    assert repr(database) == "<prolint2.QueryProteins containing 0 atoms>"

    whole = mda.Universe.empty(n_atoms=6)
    database = QueryProteins(whole)
    assert str(database) == "<prolint2.QueryProteins containing 6 atoms>"
    assert repr(database) == "<prolint2.QueryProteins containing 6 atoms>"