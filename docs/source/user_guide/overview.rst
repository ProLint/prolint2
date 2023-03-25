********
Overview
********

Improved features
=================
The major changes in **prolint2** compared to the old version that allow them to overpass the previously mentioned limitations are the migration to the `MDAnalysis`_ ecosystem, the use of a **Cython** version of a cell list algorithm greatly inspired by the neighbors grid search implemented in **GROMACS** (including the capacity to handle PBC), and the use of highly optimized data structures  (**Scipy** matrices and **Pandas** dataframes) to store the contacts.
    
    #. Routine to calculate the distance-based contacts using a cell list fixed radius neighbors search algorithm. 

    #. It reads the frames of the trajectory completely *on-the-fly*, so it does not overload memory.

    #. It takes into account the PBC in both orthorhombic and triclinic types of simulation boxes.

    #. It automatically identifies the Protein and Lipid groups for the calculation of the contacts, so you do not need to make any previous cleaning steps in your system.

    #. User-friendly, easy-to-install and well-documented tool, based on the actively-maintained `MDAnalysis`_ package.

Current performance
===================

    

.. _MDAnalysis: https://www.mdanalysis.org
