Getting Started
===============

This page details how to get started with UFCC. 

Creating the UFCC object:
-------------------------
The first step to use the **ufcc** tool is to create an *UFCC* object using a structure/topology file and a trajectory file.
All the supported formats are listed `here`_.

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('structure.gro', 'trajectory.xtc') 

By default *UFCC* will automatically identify the proteins and the membrane in the systems. For the proteins *UFCC* will identify all atoms that belong 
to a standard set of residues based in a hard-coded set of residue names (it may not work for esoteric residues). For the membrane it will identify all the lipids 
listed below:

Supported lipid types: `POPC, DPPC, DOPC, CHOL, CHL1, POPS, POPE`

But if you have other lipid types in your membrane, you can also add them at the time of creating the *UFCC* object.

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('structure.gro', 'trajectory.xtc', add_lipid_types = ['POPI', 'PSM']) 

Once you have created the *UFCC* object you will be able to access information in your query proteins and your membrane database.

.. code-block:: python

    # for query proteins
    target_system.query # ufcc.QueryProteins object

    target_system.query.whole # the whole AtomGroup from where you can select the
                              # query for the contact calculation. (this is unmutable)

    target_system.query.AG # the query AtomGroup that you will use for the contacts calculation 
                           # (it is a subselection of the target_system.query.whole AtomGroup)
                           # and can be modified anytime using target_system.query.select()

    target_system.query.list_proteins() # returns a numpy array with the different proteins in 
                                       # your target_system.query.whole AtomGroup.
                                       # The different proteins were selected based on the different
                                       # segments definid in your structure file (only work for some formats)
                                       # or otherwise assuming that proteins are ordered and the start residue 
                                       # of the next protein is always smaller than the last residue of 
                                       # the previous protein. These labels will be added in form of a metadata
                                       # to each residue and you will be able to use them for your selections.

    # for database lipids 
    target_system.database # ufcc.MembraneDatabase object

    target_system.database.whole # the whole AtomGroup from where you can select the
                              # databse for the contact calculation. (this is unmutable)

    target_system.database.AG # the database AtomGroup that you will use for the contacts calculation 
                           # (it is a subselection of the target_system.database.whole AtomGroup)
                           # and can be modified anytime using target_system.database.select()

    target_system.database.lipid_types() # returns a numpy array with the 
                                         # names of all lipids that will be analyzed.

    target_system.database.lipid_count() # returns a dictionary with the name and count of 
                                         # each lipid that will be analyzed. 

Notice that for both **query** and **database** you have two groups **whole** and **AG**. These groups
are of type MDAnalysis.core.groups.AtomGroup, so you can access from there all the topological information
that these kind of groups can offer you, including atom names, residue names, indices, the metadata included
that was defined as *macros*, etc.

Selecting the **query** and the **database**:
---------------------------------------------
To select the atomgroups for the contacts calculation you can use the **select()** method in both
**QueryProteins** and **MembraneDatabase** objects. The selection parameter can be any of:

#. an MDAnalysis Atom, Residue or AtomGroup. 
#. a string selection using the MDAnalysis selection syntax.
#. a mask using the *macros* metadate added previously. 

The last option above is very useful for selecting individual proteins as the query for the contact calculation, 
as you can use any of the labels in target_system.query.list_proteins().

.. code-block:: python

    selection_mask = target_system.query.whole.macros == 'protein0'
    target_system.query.select(selection_mask)

Getting the contacts:
---------------------
All the information of the contacts between the **query** and the **database** will be managed using the 
**target_system.contacts** object of the **Contacts** class that is automatically initializated at the beguinning.

.. code-block:: python

    target_system.contacts # ufcc.Contacts object

    target_system.contacts.contacts # None if you have not computed or loaded any contact.
                                    # Otherwise it is a numpy array of scipy.sparse matrices.

Previous to the computation of the contacts you can define the backend that you prefer using the 
the **runner** attribute of the **Contacts** class, which is an instance of the **Runner** class.
For now the **Runner** class has only two attributes *backend* and *n_jobs*, but the idea is to make 
it more complex to be able to configure the *distributed* scheduler of **Dask** to run parallel jobs 
on remote machines.

.. code-block:: python

    target_system.contacts.runner.backend # 'serial' or 'parallel'. ('serial by default')

    target_system.contacts.runner.n_jobs # number of CPU cores to use. (-1 by default, all CPU cores)

To compute the contacts you can use the **compute()** method defining the distance cutoff (in Angstroms) that you want to use 
for the contacts determination (by default 7 Angstroms).

.. code-block:: python

    target_system.contacts.compute(cutoff=7) # this will populate target_system.contacts.contacts

Save/load contacts:
-------------------
You can save/load contacts information using the **save()** and **load()** methods as below:

.. code-block:: python

    target_system.contacts.save('contacts.pkl') # this will save a pkl file with the contacts information 
                                                # stored in target_system.contacts.contacts (useful when 
                                                # you want to use the contacts information for later).

    target_system.contacts.load('old_contacts.pkl') # this will load the contacts information in a pkl file  
                                                # to target_system.contacts.contacts (useful when you have
                                                # precomputed contacts information).


Counting contacts:
-------------------
To count the contacts from the **numpy array of scipy.sparse matrices** stored in the *contacts* attribute
you can use the **count_contacts()** method, which populates the *counts* attribute.

.. code-block:: python

    target_system.count_contacts() # populates the target_system.contacts.counts attribute

    target_system.counts # None if you have not used the count_contacts() method.
                         # Otherwise it is a pandas DataFrame with the counted contacts.

.. _`here`: https://userguide.mdanalysis.org/stable/formats/index.html