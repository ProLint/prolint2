=========================================================
Welcome to the ProLint2's documentation!
=========================================================

**Last updated:** |today|

.. start-badges

|pypi|
|license|
|actions|
|codecov|
|docs|

.. |pypi| image:: https://img.shields.io/pypi/v/prolint2?color=blue
     :alt: PyPI
     :target: https://pypi.org/project/prolint2/

.. |license| image:: https://img.shields.io/badge/license-MIT-blue.svg
    :alt: license
    :target: https://opensource.org/

.. |actions| image:: https://github.com/ProLint/prolint2/workflows/CI/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/ProLint/prolint2/actions?query=workflow%3ACI

.. |codecov| image:: https://codecov.io/gh/ProLint/prolint2/graph/badge.svg
    :alt: codecov
    :target: https://app.codecov.io/gh/ProLint/prolint2

.. |docs| image:: https://readthedocs.org/projects/prolint2/badge/?version=latest
    :target: https://prolint2.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
          
.. end-badges

ProLint2 calculates distance-based lipid-protein interactions from molecular dynamics trajectories of membrane protein systems. 

.. ..  figure:: _static/logo.png
..     :align: center

Installation
============
To install **prolint2** we recommend creating a new conda environment as follows:

.. code-block:: bash

   conda create -n prolint2 python=3.8
   conda activate prolint2

Then you can install **prolint2** via pip:

.. code-block:: bash

   pip install prolint2

Basic examples:
===============

Using the Prolint2's API:

.. code-block:: python

      from prolint2 import PL2
      from prolint2.sampledata import GIRK

      target_system = PL2(GIRK.coordinates, GIRK.trajectory) 

      target_system.contacts.compute(cutoff=7)
      target_system.contacts.export('results.csv')

      
Using the Prolint2's command-line interface:

.. code-block:: bash

      prolint2 coordinates.gro trajectory.xtc -c 7 -e results.csv

You can find more details on how to use **prolint2** in the `usage page`_.

.. Contents
.. ========

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   source/intro.rst
   source/usage.rst
   api
   source/contributing.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


License 
=======

Source code included in this project is available under the `MIT License`_.

Copyright
=========
Copyright (c) 2022, Daniel P. Ramirez & Besian I. Sejdiu


Acknowledgements
================ 
The respository structure of **ProLint2** is based on the `Computational Molecular Science Python Cookiecutter <https://github.com/molssi/cookiecutter-cms>`__ version 1.6.


.. _`MIT License`: https://opensource.org/licenses/MIT
.. _`github.com/Prolint/prolint2`: https://github.com/ProLint/prolint2
.. _`usage page`: source/usage.html
