Connecting with Prolint (temporal)
==================================
As the UFCC package is still in a very early stage of development, we have mainly focused on optimizing the speed and the 
memory usage of the distance-based contacts routines, as well as creating a suitable platform to allow for the incremental
development of the package. That is why **ufcc** does not cover yet the visualization side for the analysis of the calculated 
contacts. To fulfill this, we created an **export_to_prolint()** method in the **Contacts** class that will allow for 
the visualization of the contacts using the prolint's tools. This method converts the contacts stored in the 
**ufcc**'s data structures for those used in Prolint. This is a temporary solution to those users who want to take
advantage of the fast calculation of the contacts in **ufcc**, and visualize the results in the way that they are used to in
the Prolint server or using the **prolintpy** library.

.. warning::
    The link to the Prolint's tools is nothing but natural, as **ufcc** is born inside the Prolint ecosystem, but it is important to mention
    that **ufcc** is aimed to be a much more optimized tool than what the prolint's tools currently are, that is the reason why the **export_to_prolint()**
    method can take a considerable time compared with the speed of the rest of the methods in **ufcc**. Even though, the use of this temporary
    solution (exporting the **ufcc** contacts information to the analysis tools in Prolint) is several times faster than using ProLint alone. 

The tutorial below is a brief introduction on how to use **ufcc** in combination with **prolintpy**, to get more information on how to use the **prolintpy** library
you can refer to the `official documentation`_ of this tool.

Installing **ufcc** and **prolintpy** on the same environment:
--------------------------------------------------------------

.. code-block:: bash

      conda create -n ufcc_prolint python=3.7 #(or higher)
      conda activate ufcc_prolint

      conda install -c conda-forge mdtraj
      pip install prolintpy

      pip install ufcc

.. .. warning::
..     This installation to use both packages at the same time has two issues from the **prolintpy** side:

..     * `Import error #6`_. 
..     * Change **isinstance** to **hasattr** in the **retrieve_contacts()** function of the **compute_contacts.py** file in **prolintpy**.

1. Getting the contacts information with **ufcc**:
--------------------------------------------------

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('coordinates.gro', 'trajectory.xtc') 

      target_system.contacts.compute(cutoff=7)
      target_system.contacts.count_contacts()


2. Exporting the results in prolint format to a pkl file:
---------------------------------------------------------

.. code-block:: python

      target_system.contacts.export_to_prolint(path='prolint_results.pkl')


3. Configuring system with **prolintpy**:
-----------------------------------------

.. code-block:: python

      import mdtraj as md
      import prolintpy as pl

      t = md.load('trajectory.xtc', top='coordinates.gro')

      # Specify the resolution of the input data: martini or atomistic
      resolution = "martini"

      # Define the lipid topology
      lipids = pl.Lipids(t.topology, resolution=resolution)

      # Define the protein topology
      proteins = pl.Proteins(t.topology, resolution=resolution).system_proteins()


4. Importing the contacts results to use with prolintpy:
--------------------------------------------------------

.. code-block:: python

      import pickle
    
      # loading the file created during the step 1.
      with open('prolint_results.pkl', 'rb') as f:
          results = pickle.load(f)


5. Using prolintpy tools to analyze the contacts results
--------------------------------------------------------
One of the helper functions provided by **prolintpy** is contacts_dataframe which builds a pandas DataFrame for all contacts. 
This is useful, since many of the visualization applications rely on this dataframe structure. Using this function is 
straightforward once you have the contacts results.

.. code-block:: python

      df = pl.contacts_dataframe(results, proteins, t, radius=0.7, resolution='martini')

Then you will be able to use the visualization tools of **prolintpy** as explained `here`_.

.. _`official documentation`: https://prolint.github.io/prolintpy/#/
.. _`here`: https://prolint.github.io/prolintpy/#/visualization
.. _`Import error #6`: https://github.com/ProLint/prolintpy/issues/6 