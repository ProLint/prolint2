Tutorial #5: Building Complex Contact Schemas 
==============================================

In this tutorial, we'll explore how to build complex contact representations in **prolint2**, allowing you to create intricate contact schemas for in-depth analysis of lipid-protein interactions in your molecular dynamics simulations. 

Prerequisites
-------------

Ensure you have **prolint2** installed in your Python environment.

.. code-block:: python

    from prolint2 import Universe
    from prolint2.sampledata import GIRKDataSample

    GIRK = GIRKDataSample()
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    u.normalize_by = 'actual_time'  # Choose your preferred normalization method

Overview of Complex Contact Schemas
-----------------------------------

**prolint2**'s modular design allows you to create complex contact schemas by combining contact representations computed at different cutoff values. We will demonstrate how to achieve this, showcasing its power and flexibility.

Compute Contacts at Different Cutoffs
-------------------------------------

Let's start by computing contact representations at three different cutoff values, `c1`, `c2`, and `c3`.

.. code-block:: python

    c1 = u.compute_contacts(cutoff=6)
    c2 = u.compute_contacts(cutoff=7)
    c3 = u.compute_contacts(cutoff=8)

These contact representations will serve as our building blocks for constructing the contact schemas.

Using Double Cutoffs
--------------------

We can use a double cutoff as a damping layer to measure contact durations more accurately. This is done by taking the intersection of two contact representations, `c1` and `c2`.

.. code-block:: python
    
    c3 = c1.intersection(c2)

The result, `c3`, represents the contacts that exist in both `c1` and `c2`. This intuitive approach allows you to build complex contact schemas effortlessly.

Understanding the Logic
-----------------------

The logic behind this approach is straightforward. When dealing with radial cutoffs, you can think of it as a subset relationship. In other words, if `c1` represents contacts with a smaller cutoff, it is a subset of `c2`, which represents contacts with a larger cutoff. Here's a quick summary of set operations based on this relationship:

* Union (:math:`c1 \cup c2`): All elements of `c1` are already present in `c2`, so the union results in `c2`.

* Intersection (:math:`c1 \cap c2`): All elements of `c1` are present in `c2`, so the intersection includes all elements of `c1`.

* Difference (:math:`c2 - c1`): This is the set of elements that are in `c2` but not in `c1`, representing elements unique to `c2`.

* Difference (:math:`c1 - c`): This results in the empty set (:math:`\emptyset`), as all elements of `c1` are present in `c2`.

* Symmetric Difference (:math:`c1 \bigtriangleup c2`): This is equivalent to the difference (:math:`c2 - c1`), as it includes elements unique to `c2`.

Constructing More Complex Schemas
---------------------------------

You can use the same approach to create contact schemas with more cutoffs. For example, to analyze how contacts with different lipids change with the cutoff, you can combine multiple cutoffs.

.. code-block:: python

    _ = c1 + c2  # Double cutoff
    _ = c1 + c2 + c3  # Triple cutoff
    _ = c1 + c2 + c3 + c4  # Quadruple cutoff
    _ = c1 + c2 + c3 + c4 + c5  # Quintuple cutoff

This way, you can explore how the contacts with different lipids evolve as you adjust the cutoff value.

Annular Shell Contacts
----------------------

Conversely, you can construct contact schemas that resemble a donut shape. This is useful for studying annular lipids, a critical concept in lipid-protein interactions. To do this, you subtract the core contacts from the shell contacts.

.. code-block:: python
    
    annular_shell = c2 - c1  # or c2.difference(c1)

This analysis allows you to investigate the annular lipids surrounding the protein as a function of the cutoff. You can apply this logic to multiple shells to create more intricate schemas.

.. code-block:: python

    first_shell = c2 - c1  # or c2.difference(c1)
    second_shell = c3 - c2  # or c3.difference(c2)
    third_shell = c4 - c3  # or c4.difference(c3)
    fourth_shell = c5 - c4  # or c5.difference(c4)

By combining and subtracting different contact representations, you can create complex contact schemas tailored to your specific analysis needs.
