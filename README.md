ProLint v2: an optimized tool for the analysis of lipid protein interactions.
=============================================================================

[//]: # (Badges)  
[![PyPI](https://img.shields.io/pypi/v/prolint2?color=blue)](https://pypi.org/project/prolint2/)
[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/)
[![GitHub Actions Build Status](https://github.com/ProLint/prolint2/workflows/CI/badge.svg)](https://github.com/ProLint/prolint2/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/ProLint/prolint2/graph/badge.svg)](https://app.codecov.io/gh/ProLint/prolint2)
[![docs](https://readthedocs.org/projects/prolint2/badge/?version=latest)](https://prolint2.readthedocs.io/en/latest/?badge=latest)

ProLint2 calculates the distance-based contacts between two references from molecular dynamics simulations. 

Installation
============
To install **prolint2** we recommend creating a new conda environment as follows:

``` bash
   conda create -n prolint2 python=3.8 
   conda activate prolint2
```

Then you can install **prolint2** via pip:

``` bash
   pip install prolint2
```

Basic example (from the command-line):
======================================
``` bash 
   prolint2 coordinates.gro trajectory.xtc -c 7 
```

The previous line calculates the contacts using a cutoff of 7 Angstroms between all lipids and proteins in the system, and give you a link to launch the visualization dashboard on a browser.

**prolint2** selection interface: 
=============================
You can enter in the interactive selection mode by selecting the `-i` option.

**prolint2** operates on atom selection called groups. **prolint2** defines two refence groups: the *Query* and the *Database*.
**prolint2** calculates contacts between the atoms in the *Query* group and the atoms in the *Database* group at the level of residue. By default, **prolint2** will extract proteins into the *Query* groups, and lipids into the *Database* group (if your membrane has lipids with resnames different to those recognized by default by **prolint2** [POPC, DPPC, DOPC, CHOL, CHL1, POPS, POPE] you can add them using the `-al` option). The current interface allows you to customize these selections based on your needs. We have created a few default selections and user-friendly actions that you can perform on them to reach a high level of especificity in the selections. Note that when you have created a group with the especific selection that you need, you can update either the *Database* or *Query* groups, as follow:
``` 
db Group_ID
```
or
```
qr Group_ID
```

You can see a detailed  list of the action keys below:
``` 
db  :  update database group for contacts calculation. (i.e. >> db Group_ID) where Group_ID is the identifier of the group on the list of groups.
qr  :  update query group for contacts calculation. (i.e. >> qr Group_ID) where Group_ID is the identifier of the group on the list of groups.
gb :  create subgroups grouped by any topology attribute ("names", "resnames", "masses", etc). (i.e. >> gp Group_ID resnames) creates new groups for every different resname in the Group_ID selected.
sl  :  create new groups by splitting a previous group at the desired level ("segment", "residue", "atom", "molecule"). (i.e. >> sl Group_ID residue) splits the Group_ID selected into all the residues that make it up, and create a new group for each of them.
add:  merge two or more groups and create a new group with the combination. (i.e. >> add Group_ID1 Group_ID2 ... Group_IDn) creates a new group for the combination of the groups Group_ID1 Group_ID2 ... Group_IDn.
del:  delete a group from the list of groups. (i.e. >> del Group_ID) remove the Group_ID selected.
lg :  print the list of groups.
h  :  print the help for all the available action keys.
e  :  exit interactive selection and calculate de contacts between the Query and Database groups. 
```
As for explaining the capabilities of this selection interface we will show some examples below using a membrane made up POPS, CHOL and POPE and a single protein:

**Example 1:** Using only two types of protein residues (ARG and LEU) as *Query* and all the lipids in the membrane as *Database*:
- Default groups:
```
Database and Query groups:

Database --> 7824 atoms
Query --> 5216 atoms
-----------------------------------------------------

Selection groups:

(0) : System --> 23820 atoms
(1) : Protein --> 2956 atoms
(2) : Lipids --> 20864 atoms
(3) : POPS --> 7824 atoms
(4) : CHOL --> 5216 atoms
(5) : POPE --> 7824 atoms
```
- Grouping by resnames the Protein group:
```
gb 1 resnames
```
```
(0) : System --> 23820 atoms
(1) : Protein --> 2956 atoms
(2) : Lipids --> 20864 atoms
(3) : POPS --> 7824 atoms
(4) : CHOL --> 5216 atoms
(5) : POPE --> 7824 atoms
(6) : ALA grouped by resnames from Protein --> 52 atoms
(7) : GLU grouped by resnames from Protein --> 224 atoms
(8) : THR grouped by resnames from Protein --> 192 atoms
(9) : PHE grouped by resnames from Protein --> 352 atoms
(10) : ASP grouped by resnames from Protein --> 112 atoms
(11) : HIS grouped by resnames from Protein --> 112 atoms
(12) : GLN grouped by resnames from Protein --> 80 atoms
(13) : TYR grouped by resnames from Protein --> 160 atoms
(14) : ARG grouped by resnames from Protein --> 192 atoms
(15) : ASN grouped by resnames from Protein --> 104 atoms
(16) : SER grouped by resnames from Protein --> 144 atoms
(17) : MET grouped by resnames from Protein --> 96 atoms
(18) : GLY grouped by resnames from Protein --> 76 atoms
(19) : LYS grouped by resnames from Protein --> 168 atoms
(20) : TRP grouped by resnames from Protein --> 140 atoms
(21) : VAL grouped by resnames from Protein --> 208 atoms
(22) : LEU grouped by resnames from Protein --> 240 atoms
(23) : ILE grouped by resnames from Protein --> 176 atoms
(24) : CYS grouped by resnames from Protein --> 64 atoms
(25) : PRO grouped by resnames from Protein --> 64 atoms
```
- Merging ARG and LEU groups:
```
add 14 22
```
```
(0) : System --> 23820 atoms
(1) : Protein --> 2956 atoms
(2) : Lipids --> 20864 atoms
(3) : POPS --> 7824 atoms
(4) : CHOL --> 5216 atoms
(5) : POPE --> 7824 atoms
(6) : ALA grouped by resnames from Protein --> 52 atoms
(7) : GLU grouped by resnames from Protein --> 224 atoms
(8) : THR grouped by resnames from Protein --> 192 atoms
(9) : PHE grouped by resnames from Protein --> 352 atoms
(10) : ASP grouped by resnames from Protein --> 112 atoms
(11) : HIS grouped by resnames from Protein --> 112 atoms
(12) : GLN grouped by resnames from Protein --> 80 atoms
(13) : TYR grouped by resnames from Protein --> 160 atoms
(14) : ARG grouped by resnames from Protein --> 192 atoms
(15) : ASN grouped by resnames from Protein --> 104 atoms
(16) : SER grouped by resnames from Protein --> 144 atoms
(17) : MET grouped by resnames from Protein --> 96 atoms
(18) : GLY grouped by resnames from Protein --> 76 atoms
(19) : LYS grouped by resnames from Protein --> 168 atoms
(20) : TRP grouped by resnames from Protein --> 140 atoms
(21) : VAL grouped by resnames from Protein --> 208 atoms
(22) : LEU grouped by resnames from Protein --> 240 atoms
(23) : ILE grouped by resnames from Protein --> 176 atoms
(24) : CYS grouped by resnames from Protein --> 64 atoms
(25) : PRO grouped by resnames from Protein --> 64 atoms
(26) : Group combined from ['14', '22'] --> 432 atoms
```
- Updating the *Query* with the previously created group:
```
qr 26
```
- Exiting the interactive selection mode:
```
e
```

**Example 2:** Using only a single bead (PO4) of a single type of lipid in the membrane (POPS) as *Database*, and the whole protein as *Query*:
- Default groups:
```
Database and Query groups:

Database --> 7824 atoms
Query --> 5216 atoms
-----------------------------------------------------

Selection groups:

(0) : System --> 23820 atoms
(1) : Protein --> 2956 atoms
(2) : Lipids --> 20864 atoms
(3) : POPS --> 7824 atoms
(4) : CHOL --> 5216 atoms
(5) : POPE --> 7824 atoms
```
- Grouping by names the POPS group:
```
gb 1 names
```
```
(0) : System --> 23820 atoms
(1) : Protein --> 2956 atoms
(2) : Lipids --> 20864 atoms
(3) : POPE --> 7824 atoms
(4) : POPS --> 7824 atoms
(5) : CHOL --> 5216 atoms
(6) : D2A grouped by names from POPS --> 652 atoms
(7) : C4B grouped by names from POPS --> 652 atoms
(8) : C1B grouped by names from POPS --> 652 atoms
(9) : C1A grouped by names from POPS --> 652 atoms
(10) : C3A grouped by names from POPS --> 652 atoms
(11) : GL1 grouped by names from POPS --> 652 atoms
(12) : PO4 grouped by names from POPS --> 652 atoms
(13) : C4A grouped by names from POPS --> 652 atoms
(14) : GL2 grouped by names from POPS --> 652 atoms
(15) : CNO grouped by names from POPS --> 652 atoms
(16) : C3B grouped by names from POPS --> 652 atoms
(17) : C2B grouped by names from POPS --> 652 atoms
```
- Updating the *Database* with the group corresponding the bead PO4:
```
db 12
```
- Exiting the interactive selection mode:
```
e
```

How to contribute?
==================
If you find a bug in the source code, you can help us by submitting an issue to our [GitHub repo](https://github.com/ProLint/prolint2/tree/demo). Even better, you can submit a Pull Request with a fix. 

We really appreciate your feedback!

License 
=======

Source code included in this project is available under the [MIT License](https://opensource.org/licenses/MIT).

Copyright
=========
Copyright (c) 2022, Daniel P. Ramirez & Besian I. Sejdiu


Acknowledgements
================ 
The respository structure of **ProLint2** is based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
