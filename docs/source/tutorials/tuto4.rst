Tutorial #4: Residence Time Calculations
========================================

In this tutorial, we will explore how to compute residence times for specific interactions. Residence time calculations are crucial for understanding the stability of molecular interactions in your simulations.

Prerequisites
-------------

Before we begin, make sure you have **prolint2** installed in your Python environment.

.. code-block:: python
    from prolint2 import Universe
    from prolint2.metrics.restime import KoffCalculator
    from prolint2.sampledata import GIRKDataSample

    GIRK = GIRKDataSample()
    u = Universe(GIRK.coordinates, GIRK.trajectory)
    u.normalize_by = 'actual_time'  # Ensure we use the true time for normalization
    c = u.compute_contacts(cutoff=7)

Understanding Residence Time Calculations
-----------------------------------------

Residence time calculations involve determining how long a specific molecular interaction persists in your simulation. In this tutorial, we will calculate the residence time for residue 401 in contact with cholesterol (CHOL).

Setting Up Parameters
---------------------

To compute residence time, we need to define some parameters related to the simulation data. These parameters include the total simulation time and the time step.

.. code-block:: python

    # Calculate total simulation time and time step
    totaltime = u.trajectory.totaltime * u.params['unit_conversion_factor']
    timestep = round(u.trajectory.dt * u.params['unit_conversion_factor'], 4)

    print("Total Simulation Time: {} units".format(totaltime))
    print("Time Step: {} units".format(timestep))

These parameters are essential for the residence time calculation.

Computing Residence Time
------------------------

Now, let's compute the residence time for residue 401 in contact with cholesterol (CHOL). We will use the `KoffCalculator` class from the `prolint2.metrics.restime` module.

.. code-block:: python

    # Get the data for residue 401 with CHOL
    data = c.contacts[401]['CHOL']

    # Compute residence time using KoffCalculator
    r401_chol = KoffCalculator(data, totaltime, timestep, fitting_func_name='bi_expo')

    print("Residence Time (koff): {} units".format(r401_chol.koff))
    print("Residence Time (tau): {} units".format(r401_chol.res_time))

In this example, we use the `bi_expo` fitting function to calculate residence time. The `KoffCalculator` class takes the contact data, total simulation time, time step, and fitting function as inputs to compute the residence time.
