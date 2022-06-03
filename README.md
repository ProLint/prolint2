Ultra-Fast Contacts Calculation (UFCC)
=====================================

[//]: # (Badges)
[![PyPI](https://img.shields.io/pypi/v/ufcc?color=blue)](https://pypi.org/project/ufcc/)
[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/)
[![GitHub Actions Build Status](https://github.com/ProLint/ufcc/workflows/CI/badge.svg)](https://github.com/ProLint/ufcc/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/ProLint/ufcc/graph/badge.svg)](https://app.codecov.io/gh/ProLint/ufcc)
[![docs](https://readthedocs.org/projects/ufcc/badge/?version=latest)](https://ufcc.readthedocs.io/en/latest/)




TThe **Ultra-Fast Contact Calculation (UFCC)** tool calculates the distance-based contacts between two references from molecular dynamics simulations. This release of **ufcc** is done as a concept test covering only the analysis of lipid-protein interactions on the framework 
of the Canadian Chemistry Conference and Exhibition 2022, but is aimed also to include other types of interactions in the future (i.e. protein-protein, protein-ligand, and membrane-ligand interactions).

Installation
============
To install **ufcc** we recommend creating a new conda environment as follows:

``` bash
   conda create -n ufcc python=3.7 #(or higher)
   conda activate ufcc
```

Then you can install **ufcc** via pip:

``` bash
   pip install ufcc
```


Basic example:
==============

For serial contacts calculation:

``` python
   from ufcc import UFCC

   target_system = UFCC('coordinates.gro', 'trajectory.xtc') 

   target_system.contacts.compute()
   target_system.contacts.count_contacts()
   target_system.contacts.counts
```
      
For parallel contacts calculation:

``` python
   from ufcc import UFCC

   target_system = UFCC('coordinates.gro', 'trajectory.xtc') 
   target_system.contacts.runner.backend = 'parallel'
   
   target_system.contacts.compute()
   target_system.contacts.count_contacts()
   target_system.contacts.counts
```

You can find more details on how to use **ufcc** in the [documentation page](https://ufcc.readthedocs.io/en/latest/index.html).

License 
=======

Source code included in this project is available under the [MIT License](https://opensource.org/licenses/MIT).

Copyright
=========
Copyright (c) 2022, Daniel P. Ramirez & Besian I. Sejdiu


Acknowledgements
================ 
The respository structure of **UFCC** is based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.