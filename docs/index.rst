=========================================================
Welcome to the UFCC's documentation!
=========================================================

**Last updated:** |today|

.. start-badges

|mdanalysis|
|actions|
|codecov|
|docs|

.. |mdanalysis| image:: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
    :alt: Powered by MDAnalysis
    :target: https://www.mdanalysis.org

.. |actions| image:: https://github.com/ProLint/ufcc/workflows/CI/badge.svg
    :alt: GitHub Actions Build Status
    :target: https://github.com/ProLint/ufcc/actions?query=workflow%3ACI

.. |codecov| image:: https://codecov.io/gh/ProLint/UFCC/branch/master/graph/badge.svg
    :alt: codecov
    :target: https://codecov.io/gh/ProLint/UFCC/branch/master

.. |docs| image:: https://readthedocs.org/projects/ufcc/badge/?version=latest
          :alt: docs
          :target: https://ufcc.readthedocs.io/en/latest/
          
.. end-badges

The **Ultra-Fast Contact Calculation (UFCC)** tool calculates the distance-based contacts between two references using as input the trajectories of molecular dynamics simulations. This release of **ufcc** is done as a concept test covering only the analysis of lipid-protein interactions on the framework 
of the Canadian Chemistry Conference and Exhibition 2022, but it is thinked to also include other types of interactions in the future (i.e. protein-protein, ligand-protein and ligand-membrane interactions).

.. ..  figure:: _static/logo.png
..     :align: center

Installation
============
To install **ufcc** we highly recommend to create a new conda environment as follows:

.. code-block:: bash

   conda create -n ufcc 
   conda activate ufcc

Then you can install **ufcc** using the conda-forge channel:

.. code-block:: bash

   conda install -c conda-forge ufcc

or via pip:

.. code-block:: bash

   pip install ufcc

Basic example:
==============

For serial contacts calculation:

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('coordinates.gro', 'trajectory.xtc') 

      target_system.contacts.compute()
      target_system.contacts.count_contacts()

      
For parallel contacts calculation:

.. code-block:: python

      from ufcc import UFCC

      target_system = UFCC('coordinates.gro', 'trajectory.xtc') 
      target.contacts.runner.backend = 'parallel'
      
      target_system.contacts.compute()
      target_system.contacts.count_contacts()

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
The respository structure of **UFCC** is based on the `Computational Molecular Science Python Cookiecutter <https://github.com/molssi/cookiecutter-cms>`__ version 1.6.


.. _`MIT License`: https://opensource.org/licenses/MIT
.. _`github.com/Prolint/ufcc`: https://github.com/ProLint/ufcc
.. _`usage page`: source/usage.html
