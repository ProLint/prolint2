********
Overview
********

Welcome to **prolint2**, a powerful tool for analyzing and visualizing lipid-protein interactions that builds upon its predecessor. **prolint2** comes with several new features and improved performance that enable users to overcome previous limitations. The significant enhancements incorporated in this updated version are primarily centred around a migration to the `MDAnalysis`_ ecosystem, the use of more efficient algorithms to calculate contacts based on a cutoff distance, and the implementation of new types of analysis and highly interactive visualizations to explore the resulting interactions.

Key features of **prolint2** include:

1. **Advanced Contact Calculation**: **prolint2** offers a robust routine for computing distance-based contacts using a cell list algorithm. This approach significantly enhances the scalability and efficiency of contact calculations, which becomes relevant when dealing with complex membrane systems.

2. **On-the-Fly Frame Processing**: One of the standout features of **prolint2** is its ability to read trajectory frames on-the-fly. This approach minimizes memory usage and ensures that the software remains performant even when handling extensive trajectories.

3. **Periodic Boundary Conditions (PBC) Handling**: **Prolint2** is able to handle both orthorhombic and triclinic types of simulation boxes, taking into account periodic boundary conditions (PBC) to provide accurate results in various molecular dynamics simulations.

4. **Automated Group Identification**: This version of **prolint2** boasts automatic identification of protein and lipid groups within the system. Users no longer need to perform extensive manual cleaning and grouping before initiating contact calculations, streamlining the workflow.

5. **High Modularity for Custom Analysis**: **Prolint2** is characterized by a high level of modularity, enabling users to define new types of analyses and interaction metrics. Researchers can tailor their analyses to specific research questions, facilitating in-depth investigations of lipid-protein interactions in a manner that suits their unique needs and objectives.

6. **Arithmetic Functions for Contact Groups**: The software includes built-in arithmetic functions that allow users to perform basic calculations between contact groups. This feature is especially valuable for studying annular distributions of lipids around proteins. 

7. **Interactive Multi-Application Dashboard**: **Prolint2** introduces a highly interactive multi-application dashboard that offers an intuitive and user-friendly interface for exploring the results of lipid-protein interactions. This dashboard streamlines the visualization of data and facilitates the interpretation of complex interaction patterns.

8. **Residence Time Calculation**: **Prolint2** provides the ability to calculate residence times for specific lipid types interacting with particular protein residues. This feature is vital for understanding the dynamics and stability of these interactions over time, shedding light on their functional relevance.

9. **Multilevel Interaction Analysis**: Researchers can conduct multilevel analyses of interactions, including residue-lipid, residue-lipid type, residue type-lipid, and atom-atom interactions. This granular approach enables a detailed examination of interactions at various hierarchical levels, enhancing the understanding of molecular events.

10. **Binding Site Identification**: **Prolint2** includes a specialized function for identifying binding sites through shared contact analysis. This is a crucial feature for pinpointing regions where lipids and proteins form stable and significant interactions.

.. _MDAnalysis: https://www.mdanalysis.org
