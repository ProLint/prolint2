ProLint v2: an optimized tool for the analysis of lipid protein interactions.
=============================================================================

[//]: # (Badges)  
[![PyPI](https://img.shields.io/pypi/v/prolint2?color=blue)](https://pypi.org/project/prolint2/)
[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/)
[![GitHub Actions Build Status](https://github.com/ProLint/prolint2/workflows/CI/badge.svg)](https://github.com/ProLint/prolint2/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/ProLint/prolint2/graph/badge.svg)](https://app.codecov.io/gh/ProLint/prolint2)
[![docs](https://readthedocs.org/projects/prolint2/badge/?version=latest)](https://prolint2.readthedocs.io/en/latest/?badge=latest)

ProLint2 calculates distance-based lipid-protein interactions from molecular dynamics trajectories of membrane protein systems.

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

Basic examples:
===============

Using the Prolint2's API:

``` python
   from prolint2 import PL2
   from prolint2.sampledata import GIRK

   target_system = PL2(GIRK.coordinates, GIRK.trajectory) 

   target_system.contacts.compute(cutoff=7)
   target_system.contacts.export('results.csv')
```
      
Using the Prolint2's command-line interface:

``` bash
   prolint2 coordinates.gro trajectory.xtc -c 7 -e results.csv
```

You can find more details on how to use **prolint2** in the [documentation](https://prolint2.readthedocs.io/en/latest/index.html).

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
