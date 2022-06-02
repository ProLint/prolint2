Overview
========
The fixed-radius near neighbors is an old computer science problem that can be defined in short as finding all pairs of points within a certain distance apart as represented in the figure below.

..  figure:: ../_static/frnns.png
    :align: center

There have been proposed different approaches to tackle this problem, mainly grouped in three categories: i-brute force; ii-tree based; and iii-grid based methods, which are all useful in certain range of conditions, especially
determined by the way they scale respect to the number of points in the system. 

The brute-force approach compute the distances between all pairs of points in the dataset, which for **N** samples
in **D** dimensions, scales as :math:`\textbf{O}(\textbf{DN}^{2})`. This can useful for an small number of points, but it become very slow when you 
are working with large **N**. 

The tree-based data structures attempt to reduce the required number of distance calculations by efficiently encoding aggregate distance information 
for the sample. The basic idea is that if point A is very distant from point B, and point B is very close to point C, then we know
that points A and C are very distant, *without having to explicitly calculate their distance*. In this way, the computational cost 
of a nearest neighbors search can be reduced to :math:`\textbf{O}(\textbf{DN}\log{(\textbf{N})})` or better.
The most common tree-based models for low-dimensional datasets (:math:`D<20`) 
are the KD-trees. A KD-tree is a binary tree structure which recursively partitions 
the parameter space along the data axes, dividing it into nested orthotropic regions into which data points are filed. The construction of a 
KD-tree is very fast because partitioning is performed only along the data axes, no D-dimensional distances need to be computed. Once constructed, 
the nearest neighbors of a query point can be determined with only :math:`O(log(\textbf{N}))`` distance computations. This is a significant improvement over brute-force for large **N**, 
but it is still very slow when you are treating with very large **N**. 

In the case of the grid based methods, they divide the system into smaller subdomains called cells and distributes
every particle to these cells based on their positions. Subsequently, any distance based query first
identifies the corresponding cell position in the domain followed by distance evaluations within
the identified cell and neighboring cells only (each cell has a size equal to or slightly larger than
the cutoff radius :math:`r_c`; so each particle in a given cell interacts with only those particles in the same
or neighboring cells). Since the allocation of a particle to a cell is an operation that scales with **N**
and the total number of cells that needs to be considered for the calculation of the interaction is
independent of the system (always 9 cells in a bidimensional space, or 27 for a tridimensional space), the cell list method
allow the computation of neighbors with a cost of :math:`O(\textbf{N})`, instead of :math:`O(\textbf{N}^2)`.

..  figure:: ../_static/neighbours_list.png
    :align: center

The current neighbors search implementation of **prolintpy** involves the brute-force computation of distances between 
all pairs of points in the dataset, which is very inefficient for lipid-protein systems with very high number of atoms. 
Furthermore, the current implementation only looks for the neighbors based in raw distances, whitout taking into account 
the periodic boundary conditions (PBC) present in most molecular dynamic simulations. Other limitations in the current 
version of **prolintpy** include: i- high memory usage as it needs to load all the trajectory at once to run the calculation of the contacts;
ii- need to modify the trajectories previous to the setup of the Lipid and Protein groups for the calculation of the contacts; and
iii- core dependencies no longer mantained for the handling of the molecular dynamic trajectories (MDTraj).

How **ufcc** is able to solve the limitations in **prolintpy**?
---------------------------------------------------------------

The major changes in **ufcc** compared with **prolintpy** that allows it to overpass the previously mentioned limitations is the migration
to the `MDAnalysis`_ ecosystem, the use of the a **Cython** version of a cell list algorithm greatly inspired by the
neighbors grid search implemented in **GROMACS** (including the capacity to handle PBC), and the use of highly optimized data structures 
(**Scipy** matrices and **Pandas** dataframes) to handle the contacts.

Features
--------

#. Serial and parallel routines to calculate the distance-based contacts using a cell list fixed radius neighbors search algorithm. Both routines are able to scale 
lineally with the number of atoms, and in the case of the 'parallel' routine it uses a *split-apply-combine* approach to deal with the analysis of very large trajectories, 
so there is practically no limitations in size and time respect to the simulations that can be handled. 

#. It reads the frames of the trajectory completely *on-the-fly*, so it does not need to load anything on memory.

#. It deals with PBC in both orthorombic and triclinic types of simulation boxes.

#. It automatically identify the Protein and Lipid groups for the calculation of the contacts, so you don't need to make any previous cleaning step of your system.

#. 

.. _MDAnalysis: https://www.mdanalysis.org
