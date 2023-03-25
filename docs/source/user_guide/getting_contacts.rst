Getting the contacts
====================
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

