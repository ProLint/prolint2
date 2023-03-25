Selecting the **query** and the **database**
=============================================
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

