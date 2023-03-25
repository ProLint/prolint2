***********************
Creating the PL2 object
***********************
The first step to use the **prolint2** tool is to create a *PL2* object using a structure/topology file and a trajectory file.
All the supported formats are listed `here`_.

.. code-block:: python

      from prolint2 import PL2
      from prolint2.sampledata import GIRK

      target_system = PL2(GIRK.coordinates, GIRK.trajectory) 

Supported residue types
=======================
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

.. _`here`: https://userguide.mdanalysis.org/stable/formats/index.html
