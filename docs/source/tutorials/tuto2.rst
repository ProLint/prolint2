Tutorial #2: **prolint2** `Universe` Properties 
===============================================

In this tutorial, we will explore how to access and modify the properties of the `Universe` object to customize your analysis of lipid-protein interactions. You will learn how to adjust parameters such as units, normalization methods, and more to gain a deeper understanding of your molecular dynamics simulations.

Prerequisites
-------------

Before we begin, make sure you have **prolint2** installed in your Python environment.

.. code-block:: python

    from prolint2 import Universe
    from prolint2.sampledata import GIRKDataSample
    GIRK = GIRKDataSample()
    u = Universe(GIRK.coordinates, GIRK.trajectory)

Displaying Current `Universe` Properties
----------------------------------------

Let's start by displaying the current properties of the `Universe` object using the `params` attribute. These properties include units, normalization methods, conversion factors, and more.

.. code-block:: python

    # Display the current properties of the Universe object
    u.params

Here's a breakdown of the parameters and their meanings:

* `units`: This parameter determines the unit in which the results will be displayed.
* `normalizer`: It defines how the results are normalized, with options for 'counts' (raw frame counts), 'actual_time' (normalized by the true simulation time), and 'time_fraction' (normalized by fractions of time spent in contact).
* `unit_conversion_factor`: This factor converts raw counts to the desired unit.
* `norm_factor`: This factor is updated based on the chosen normalization type.

Computing Contacts with Default Parameters
------------------------------------------

Before modifying the parameters, let's compute contacts with the default settings. Contacts are measured as contact durations relative to the total simulation time and expressed in microseconds (`us`).

.. code-block:: python

    # Compute contacts with default parameters
    default = u.compute_contacts(cutoff=7)

    # Access contacts for a specific residue (e.g., residue 17 from 'POPE')
    residue17_default_contacts = default.contacts.get(17).get('POPE')

    # Display the first 5 contact durations
    residue17_default_contacts[:5]

Changing `Universe` Parameters and Re-computing Contacts
--------------------------------------------------------

Now, let's modify the `Universe` parameters to suit your analysis needs.

1. Change the Normalization Method

We will switch from 'time_fraction' to 'actual_time' for normalization. This will represent contact durations in relation to the actual time of the trajectory.

.. code-block:: python

    # Change the normalization method to 'actual_time'
    u.normalize_by = 'actual_time'

    # Display the updated parameters
    u.params


Notice how the `norm_factor` is updated based on the chosen normalization type.

2. Change the Units

Let's change the units from microseconds (`us`) to nanoseconds (`ns`) to measure contact durations in their true simulation time.

.. code-block:: python

    # Change the units to nanoseconds
    u.units = 'ns'

    # Display the updated parameters
    u.params

Notice how both the `unit_conversion_factor` and `norm_factor` are adjusted based on the chosen units.

Re-compute Contacts with Modified Parameters
--------------------------------------------

Now, let's re-compute the contacts with the adjusted parameters.

.. code-block:: python

    # Re-compute contacts with the modified parameters
    updated = u.compute_contacts(cutoff=7)

    # Access contacts for the same residue (e.g., residue 17 from 'POPE')
    residue17_modified_contacts = updated.contacts.get(17).get('POPE')

    # Display the first 5 contact durations
    residue17_modified_contacts[:5]

Congratulations! You've successfully modified `Universe` properties to customize the analysis of lipid-protein interactions in **prolint2**. This flexibility allows you to adapt your analysis to different conditions and represent your results in the most informative way for your research or simulation requirements.