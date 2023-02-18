Usage Documentation
===================

This page details how to get started with **prolint2** to work with lipid-protein systems. On the `Prolint's resources page`_ you can get some simple lipid-protein systems to test the tools before moving on with your own systems.

Using the API:
--------------
Creating the PL2 object:
-------------------------
The first step to use the **prolint2** tool is to create a *PL2* object using a structure/topology file and a trajectory file.
All the supported formats are listed `here`_.

.. code-block:: python

      from prolint2 import PL2
      from prolint2.sampledata import GIRK

      target_system = PL2(GIRK.coordinates, GIRK.trajectory) 

By default, **prolint2** will automatically identify the proteins and the membrane in the systems. For the proteins **prolint2** will identify all atoms that belong 
to a standard set of residues based on a hard-coded set of residue names (it may not work for esoteric residues). For the membrane it will identify all the lipids 
listed below:

Supported lipid types: `DLPC, DPPC, DOPC, DIPC, POPC, DPPE, DOPE, POPE, DPPS, DOPS, POPS, DPPG, DOPG, POPG, DPPA, DOPA, POPA, DPPI, POPI, DPP1, POP1, DPP2, POP2, PODG, LPC, PPC, OPC, DPSM, POSM, DPCE, DPGS, DPG1, DPG3, DPMG, CHOL, CHL1`.

But if you have other lipid types in your membrane, you can also add them at the time of initializing the *PL2* instance.

.. code-block:: python

      from prolint2 import PL2
      from prolint2.sampledata import GIRK

      target_system = PL2(GIRK.coordinates, GIRK.trajectory, add_lipid_types = ['TOG']) 

Once you have initialized the *PL2* instance you will be able to access information in your query proteins and your membrane database, 
including an additional label that is automatically added to each residue and that we called *macros*.

The *macros* metadata is going to be useful for the selection of the query and the database groups. If the residue is included in the MembraneDatabase group, 
then the *macro* metadata will be set as **membrane**; if the residue is included in the QueryProteins group then the *macro* metadata will be set with **Protein#**
depending on the number of segments (or chains) in the system; otherwise the *macro* metadata will be set as **other**.

.. code-block:: python

    ###################### Query proteins ########################

    target_system.query # prolint2.QueryProteins object

    target_system.query.whole 
    # This is the whole AtomGroup from where you can select the
    # query for the contact calculation. (it is immutable)

    target_system.query.selected 
    # This is the selected AtomGroup that you will use as query for the 
    # calculation of the contacts. It is a subselection of the query.whole AtomGroup
    # and can be modified anytime using target_system.query.select().

    target_system.query.list_proteins() 
    # This returns a numpy array with the different proteins in 
    # your target_system.query.whole AtomGroup.
    # The different proteins were selected based on the different
    # segments defined in your structure file (only work for some formats)
    # or otherwise assuming that proteins are ordered and the start residue 
    # of the next protein is always smaller than the last residue of 
    # the previous protein. These labels will be added in form of metadata
    # to each residue and you will be able to use them for your selections.

    ##################### Membrane Database ######################

    target_system.database # prolint2.MembraneDatabase object

    target_system.database.whole 
    # This is the whole AtomGroup from where you can select the
    # database for the contact calculation. (it is immutable)

    target_system.database.selected 
    # This is the selected AtomGroup that you will use as databasefor the 
    # calculation of the contacts. It is a subselection of the database.whole AtomGroup
    # and can be modified anytime using target_system.database.select().

    target_system.database.lipid_types() 
    # This returns a numpy array with the 
    # names of all the lipids that will be analyzed.

    target_system.database.lipid_count()
    # This returns a dictionary with the name and count of 
    # each lipid that will be analyzed. 

Notice that for both **query** and **database** you have two groups (**whole** and **selected**). These groups
are of type MDAnalysis.core.groups.AtomGroup, so you can access from there all the topological information
that this kind of group can offer you, including atom names, residue names, indices, the *macros* metadata, etc.

Selecting the **query** and the **database**:
---------------------------------------------
To select the references for the calculation of the contacts you can use the **select()** method in both
**QueryProteins** and **MembraneDatabase** objects. The selection parameter can be any of:

#. an MDAnalysis Atom, Residue or AtomGroup. 
#. a string selection using the MDAnalysis selection syntax.
#. a mask using the *macros* metadata added. 

The last option above is very useful for selecting individual proteins as the query for the contact calculation, 
as you can use any of the labels in target_system.query.list_proteins().

.. code-block:: python

    selection_mask = target_system.query.whole.macros == 'Protein0'
    target_system.query.select(selection_mask)

Getting the contacts:
---------------------
All the information of the contacts between the **query** and the **database** will be handled using the 
**target_system.contacts** instance of the **Contacts** class that is automatically initialized at the beginning.

.. code-block:: python

    target_system.contacts # prolint2.Contacts object

    target_system.contacts.contacts 
    # This is None if you have not computed or loaded any contact.
    # Otherwise it is a numpy array of scipy.sparse matrices.

To compute the contacts you can use the **compute()** method defining the distance cutoff (in Angstroms) that you want to use 
for the determination of the contact (by default 7 Angstroms).

.. code-block:: python

    target_system.contacts.compute(cutoff=7) 
    # this will populate target_system.contacts.contacts

Export contacts:
-------------------
You can export contacts information using the **export** method as below:

.. code-block:: python

    target_system.contacts.export('results.csv') 
    # this will export two csv files, one with the contacts information 
    # stored in target_system.contacts.contacts and a second one 
    # ('results_metrics.csv') with the metrics information.


Using the command-line:
-----------------------

.. code-block:: none

    prolint2 coordinates.gro trajectory.xtc -c 7 

The previous line calculates the contacts using a cutoff of 7 Angstroms between all lipids and proteins in the system, and give you a link to launch the visualization dashboard on a browser.

**prolint2** selection interface: 
---------------------------------
You can enter in the interactive selection mode by selecting the `-i` option.

**prolint2** calculates contacts between the atoms in the *Query* group and the atoms in the *Database* group at the level of residue. By default, **prolint2** will extract proteins into the *Query* groups, and lipids into the *Database* group (if your membrane has lipids with resnames different to those recognized by default by **prolint2** you can add them using the `-al` option). The current interface allows you to customize these selections based on your needs. We have created a few default selections and user-friendly actions that you can perform on them to reach a high level of especificity in the selections. Note that when you have created a group with the especific selection that you need, you can update either the *Database* or *Query* groups, as follow:

.. code-block:: none 

    db Group_ID

or

.. code-block:: none 

    qr Group_ID

You can see a detailed  list of the action keys below:

.. code-block:: none  

    db  : update database group for contacts calculation. (i.e. >> db Group_ID)
          where Group_ID is the identifier of the group on the list of groups.
    qr  : update query group for contacts calculation. (i.e. >> qr Group_ID) 
          where Group_ID is the identifier of the group on the list of groups.
    gb :  create subgroups grouped by any topology attribute ("names", "resnames", "masses", etc). 
          (i.e. >> gp Group_ID resnames) creates new groups for every different resname in the 
          Group_ID selected.
    sl  : create new groups by splitting a previous group at the desired level 
          ("segment", "residue", "atom", "molecule"). (i.e. >> sl Group_ID residue) 
          splits the Group_ID selected into all the residues that make it up, 
          and create a new group for each of them.
    add:  merge two or more groups and create a new group with the combination. 
          (i.e. >> add Group_ID1 Group_ID2 ... Group_IDn) creates a new group for the combination 
          of the groups Group_ID1 Group_ID2 ... Group_IDn.
    del:  delete a group from the list of groups. (i.e. >> del Group_ID) 
          remove the Group_ID selected.
    lg :  print the list of groups.
    h  :  print the help for all the available action keys.
    e  :  exit interactive selection and calculate de contacts between the 
          Query and Database groups. 

As for explaining the capabilities of this selection interface we will show some examples below using a membrane made up POPS, CHOL and POPE and a single protein:

**Example 1:** 
--------------
Using only two types of protein residues (ARG and LEU) as *Query* and all the lipids in the membrane as *Database*:
- Default groups:

.. code-block:: none 

    Database and Query groups:

    Database --> 7824 atoms
    Query --> 5216 atoms
    -----------------------------------------------------

    Selection groups:

    (0) : System --> 23820 atoms
    (1) : Protein --> 2956 atoms
    (2) : Lipids --> 20864 atoms
    (3) : POPS --> 7824 atoms
    (4) : CHOL --> 5216 atoms
    (5) : POPE --> 7824 atoms

- Grouping by resnames the Protein group:

.. code-block:: none 

    gb 1 resnames

.. code-block:: none 

    (0) : System --> 23820 atoms
    (1) : Protein --> 2956 atoms
    (2) : Lipids --> 20864 atoms
    (3) : POPS --> 7824 atoms
    (4) : CHOL --> 5216 atoms
    (5) : POPE --> 7824 atoms
    (6) : ALA grouped by resnames from Protein --> 52 atoms
    (7) : GLU grouped by resnames from Protein --> 224 atoms
    (8) : THR grouped by resnames from Protein --> 192 atoms
    (9) : PHE grouped by resnames from Protein --> 352 atoms
    (10) : ASP grouped by resnames from Protein --> 112 atoms
    (11) : HIS grouped by resnames from Protein --> 112 atoms
    (12) : GLN grouped by resnames from Protein --> 80 atoms
    (13) : TYR grouped by resnames from Protein --> 160 atoms
    (14) : ARG grouped by resnames from Protein --> 192 atoms
    (15) : ASN grouped by resnames from Protein --> 104 atoms
    (16) : SER grouped by resnames from Protein --> 144 atoms
    (17) : MET grouped by resnames from Protein --> 96 atoms
    (18) : GLY grouped by resnames from Protein --> 76 atoms
    (19) : LYS grouped by resnames from Protein --> 168 atoms
    (20) : TRP grouped by resnames from Protein --> 140 atoms
    (21) : VAL grouped by resnames from Protein --> 208 atoms
    (22) : LEU grouped by resnames from Protein --> 240 atoms
    (23) : ILE grouped by resnames from Protein --> 176 atoms
    (24) : CYS grouped by resnames from Protein --> 64 atoms
    (25) : PRO grouped by resnames from Protein --> 64 atoms

- Merging ARG and LEU groups:

.. code-block:: none 

    add 14 22

.. code-block:: none 

    (0) : System --> 23820 atoms
    (1) : Protein --> 2956 atoms
    (2) : Lipids --> 20864 atoms
    (3) : POPS --> 7824 atoms
    (4) : CHOL --> 5216 atoms
    (5) : POPE --> 7824 atoms
    (6) : ALA grouped by resnames from Protein --> 52 atoms
    (7) : GLU grouped by resnames from Protein --> 224 atoms
    (8) : THR grouped by resnames from Protein --> 192 atoms
    (9) : PHE grouped by resnames from Protein --> 352 atoms
    (10) : ASP grouped by resnames from Protein --> 112 atoms
    (11) : HIS grouped by resnames from Protein --> 112 atoms
    (12) : GLN grouped by resnames from Protein --> 80 atoms
    (13) : TYR grouped by resnames from Protein --> 160 atoms
    (14) : ARG grouped by resnames from Protein --> 192 atoms
    (15) : ASN grouped by resnames from Protein --> 104 atoms
    (16) : SER grouped by resnames from Protein --> 144 atoms
    (17) : MET grouped by resnames from Protein --> 96 atoms
    (18) : GLY grouped by resnames from Protein --> 76 atoms
    (19) : LYS grouped by resnames from Protein --> 168 atoms
    (20) : TRP grouped by resnames from Protein --> 140 atoms
    (21) : VAL grouped by resnames from Protein --> 208 atoms
    (22) : LEU grouped by resnames from Protein --> 240 atoms
    (23) : ILE grouped by resnames from Protein --> 176 atoms
    (24) : CYS grouped by resnames from Protein --> 64 atoms
    (25) : PRO grouped by resnames from Protein --> 64 atoms
    (26) : Group combined from ['14', '22'] --> 432 atoms

- Updating the *Query* with the previously created group:

.. code-block:: none 

    qr 26

- Exiting the interactive selection mode:

.. code-block:: none 

    e

**Example 2:** 
--------------
Using only a single bead (PO4) of a single type of lipid in the membrane (POPS) as *Database*, and the whole protein as *Query*:
- Default groups:

.. code-block:: none 

    Database and Query groups:

    Database --> 7824 atoms
    Query --> 5216 atoms
    -----------------------------------------------------

    Selection groups:

    (0) : System --> 23820 atoms
    (1) : Protein --> 2956 atoms
    (2) : Lipids --> 20864 atoms
    (3) : POPS --> 7824 atoms
    (4) : CHOL --> 5216 atoms
    (5) : POPE --> 7824 atoms

- Grouping by names the POPS group:

.. code-block:: none 

    gb 1 names

.. code-block:: none 

    (0) : System --> 23820 atoms
    (1) : Protein --> 2956 atoms
    (2) : Lipids --> 20864 atoms
    (3) : POPE --> 7824 atoms
    (4) : POPS --> 7824 atoms
    (5) : CHOL --> 5216 atoms
    (6) : D2A grouped by names from POPS --> 652 atoms
    (7) : C4B grouped by names from POPS --> 652 atoms
    (8) : C1B grouped by names from POPS --> 652 atoms
    (9) : C1A grouped by names from POPS --> 652 atoms
    (10) : C3A grouped by names from POPS --> 652 atoms
    (11) : GL1 grouped by names from POPS --> 652 atoms
    (12) : PO4 grouped by names from POPS --> 652 atoms
    (13) : C4A grouped by names from POPS --> 652 atoms
    (14) : GL2 grouped by names from POPS --> 652 atoms
    (15) : CNO grouped by names from POPS --> 652 atoms
    (16) : C3B grouped by names from POPS --> 652 atoms
    (17) : C2B grouped by names from POPS --> 652 atoms

- Updating the *Database* with the group corresponding the bead PO4:

.. code-block:: none 

    db 12

- Exiting the interactive selection mode:

.. code-block:: none 

    e

.. _`here`: https://userguide.mdanalysis.org/stable/formats/index.html
.. _`Prolint's resources page`: https://www.prolint.ca/resources/data