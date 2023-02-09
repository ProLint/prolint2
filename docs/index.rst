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

.. |actions| image:: https://github.com/ProLint/ufcc/workflows/CI/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/ProLint/ufcc/actions?query=workflow%3ACI

.. |codecov| image:: https://codecov.io/gh/ProLint/ufcc/graph/badge.svg
    :alt: codecov
    :target: https://app.codecov.io/gh/ProLint/ufcc

.. |docs| image:: https://readthedocs.org/projects/prolint2/badge/?version=latest
    :target: https://prolint2.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
          
.. end-badges

The **Ultra-Fast Contact Calculation (UFCC)** tool calculates the distance-based contacts between two references from molecular dynamics simulations. This release of **ufcc** is done as a concept test covering only the analysis of lipid-protein interactions on the framework 
of the Canadian Chemistry Conference and Exhibition 2022, but is aimed also to include other types of interactions in the future (i.e. protein-protein, protein-ligand, and membrane-ligand interactions).

.. ..  figure:: _static/logo.png
..     :align: center

Installation
============
To install **ufcc** we recommend creating a new conda environment as follows:

.. code-block:: bash

   conda create -n ufcc python=3.7 #(or higher)
   conda activate ufcc

Then you can install **ufcc** via pip:

.. code-block:: bash

   pip install ufcc

Basic examples:
===============

For serial contacts calculation:

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('coordinates.gro', 'trajectory.xtc') 

      target_system.contacts.compute()
      target_system.contacts.count_contacts()
      target_system.contacts.counts

      
For parallel contacts calculation:

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('coordinates.gro', 'trajectory.xtc') 
      target_system.contacts.runner.backend = 'parallel'
      
      target_system.contacts.compute()
      target_system.contacts.count_contacts()
      target_system.contacts.counts

You can find more details on how to use **ufcc** in the `usage page`_.

.. Contents
.. ========

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   
   source/intro.rst
   source/usage.rst
   source/connecting_with_prolint.rst
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
.. _`github.com/Prolint/ufcc`: https://github.com/ProLint/ufcc
.. _`usage page`: source/usage.html
