Tutorial #1: Working with the **prolint2** `Universe` 
=====================================================

In this tutorial, we will focus on working with the **prolint2** `Universe` object, which allows you to access and manipulate data related to lipid-protein interactions. You will learn how to load your trajectory data, and how to modify the **query** and **database** groups to analyze lipid-protein interactions.

Prerequisites
-------------

Before you begin, make sure you have **prolint2** installed in your Python environment.

.. code-block:: python

    from prolint2 import Universe
    from prolint2.sampledata import GIRKDataSample

    GIRK = GIRKDataSample()

Creating a **prolint2** `Universe` Instance
-------------------------------------------

ProLint can read data directly from an **MDAnalysis** `Universe` object, which is often used in molecular dynamics simulations. You can create a **prolint2** `Universe` object from an existing **MDAnalysis** `Universe`, and then you can further manipulate the **query** and **database** AtomGroups.

* Using MDAnalysis to Create a Universe Instance

.. code-block:: python

    # Create an MDAnalysis Universe instance
    mda_u = MDUniverse(GIRK.coordinates, GIRK.trajectory)

    # Define custom query and database AtomGroups
    mda_u_query = mda_u.select_atoms('protein and name BB')
    mda_u_db = mda_u.select_atoms('resname POPE')

    # Create a ProLint Universe instance from the MDAnalysis Universe
    u = Universe(universe=mda_u)

Alternatively, you can create a ProLint Universe instance directly from the MDAnalysis `Universe` with query and database information:

.. code-block:: python

    # Create a ProLint Universe instance with query and database
    u = Universe(universe=mda_u, query=mda_u_query, db=mda_u_db)


Accessing the **Query** and **Database** AtomGroups
---------------------------------------------------

Both the **query** and **database** AtomGroups in the **prolint2** `Universe` are wrappers around the **MDAnalysis** `AtomGroups`. This allows you to leverage all the functionality of `MDAnalysis AtomGroups` while also benefiting from ProLint-specific functions for lipid-protein interaction analysis.

* Example: Accessing AtomGroup Information

You can access various information about the **query** and **database** `AtomGroups`:

.. code-block:: python

    # Accessing unique residue names in the database
    u.database.unique_resnames

    # Get residue names for specific residue numbers
    u.database.get_resnames([2345, 2346, 3050])

    # Filter residue IDs by a specific residue name
    u.database.filter_resids_by_resname([2345, 2346, 3050], 'CHOL')

    # Get residue names for specific residue numbers and store the result in a dictionary
    u.query.get_resnames([1, 2, 3, 4, 5], out=dict)

Modifying the **Query** and **Database** AtomGroups
---------------------------------------------------

**prolint2** provides methods for modifying the **query** and **database** `AtomGroups`. These modifications can be useful for tailoring the `AtomGroups` to your specific analysis needs.

* Removing from **prolint2** `AtomGroup`

You can remove elements from an AtomGroup using the `remove` method. The method combines all input arguments into a single selection string concatenated with `or` statements. The resulting AtomGroup is a new instance, and it does not modify the original AtomGroup.

.. code-block:: python

    # Remove all residues with resname 'ARG' from the query
    s = u.query.remove(resname='ARG')

    # Remove all residues with resname 'ARG' and all residue numbers lower than 100
    s = u.query.remove(resname='ARG', resnum=[*range(100)])

    # More complex example: Remove all residues with resname 'ARG' and the residue number 1, 
    # and all atoms with the name 'BB' and the atomids 1-9
    s = u.query.remove(resname='ARG', resnum=[1], atomname=['BB'], atomids=[1, 2, 3, 4, 5, 6, 7, 8, 9])

    # To modify the original AtomGroup, use assignment
    u.query = u.query.remove(resname='ARG')

* Adding to prolint2 `AtomGroup`

To add elements back to an `AtomGroup`, you can use the `add` method. This allows you to revert any removals or extend the `AtomGroup`.

.. code-block:: python

    # Add back the residues that were removed from the query
    u.query = u.query.add(resname='ARG')

Computing Contacts
------------------

After modifying the **query** `AtomGroup`, you can compute lipid-protein contacts using the `compute_contacts` method. This method calculates interactions between the query and the database based on a specified cutoff distance.

.. code-block:: python

    # Compute contacts with a cutoff distance of 7 angstroms
    c = u.compute_contacts(cutoff=7)

    # Get the number of residues for which contacts have been computed
    len(c.contact_frames.keys())

By following this tutorial, you should now be comfortable with creating and working with the **prolint2** `Universe` object for lipid-protein interaction analysis and visualization. This tool offers a flexible and powerful way to study molecular interactions in your simulations.