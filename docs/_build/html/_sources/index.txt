.. Dabble documentation master file, created by
   sphinx-quickstart on Fri Jun 26 11:29:41 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Dabble
======

Dabble is a tool for building membrane protein systems. The ultimate goal of the
project is to create an easy to use, one stop tool for system construction and
parameterization.

Dabble includes both the main program dabble.py as well as a powerful API for 
building and parameterizing molecular dynamics systems.

To get started via `Anaconda Python <https://www.continuum.io/downloads>`_,
use::
    
    conda install -c rbetz dabble

.. figure:: _static/dabblebox.png
    :align: right
    Prepared system with lipid, water, protein, ligand, and ions
    
Features
--------

- Prepare membrane protein systems by inserting them into a membrane
- Prepare solvated proteins by adding water
- Add ions to neutralize and/or to desired concentration
- Parameterize with CHARMM or AMBER parameter sets
- Outputs files for simulation with AMBER or CHARMM 
- Automatic detection of post-translational modifications
- Modified amino acids made easy
- Ligands made easy! No more messing with atom names. 

Supported Parameter Sets
------------------------

- [X] CHARMM
- [X] AMBER
- [ ] Gromacs

Supported Simulation Programs
-----------------------------

- [X] CHARMM
- [X] AMBER
- [X] Anton (via conversion)
- [ ] Gromacs

Contributing
------------

Dabble is a written by Robin Betz, in the Ron Dror lab at Stanford University.
Bug finding is always apprecciated, as well as corner cases where your protein
won't dabble.

Dabble is licensed under the GPLv2 license.

- Issue tracker: https://github.com/drorlab/dabble/issues
- Source code: https://github.com/drorlab/dabble/project

Usage
-----

- :ref:`command_line`
- :ref:`builder_api`
- :ref:`parameter_api`
- :ref:`utility_programs`

API Documentation
-----------------

.. toctree::
   :maxdepth: 2
   :hidden:

   command_line
   builder_api
   parameter_api
   utility_programs
   Dabble.rst
   DabbleUtils.rst

* :ref:`search`

