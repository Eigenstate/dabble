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

Features
--------

- Prepare membrane protein systems by inserting them into a membrane
- Prepare solvated proteins by adding water
- Add ions to neutralize and/or to desired concentration
- Parameterize with CHARMM parameter sets
- Outputs files for simulation with AMBER or CHARMM 
- Automatic detection of post-translational modifications
- Modified amino acids made easy
- Ligands made easy  

Supported Parameter Sets
------------------------

- [X] CHARMM
- [ ] AMBER
- [ ] Gromacs

Supported Simulation Programs
-----------------------------

- [X] CHARMM
- [X] AMBER
- [X] Anton (via conversion)
- [ ] Gromacs

Installation
------------

A conda installer for Dabble is coming soon.
Dabble's main dependency is vmd-python.

    pip install -i https://pypi.anaconda/org/rbetz/simple vmd-python

Then, download the source code from our github.

    git clone https://github.com/drorlab/dabble.git

And install

    cd dabble
    python setup.py install

Usage
-----

.. toctree::
   :maxdepth: 2


   Usage.rst


Contribute
----------

Dabble is a project of the Ron Dror Lab at Stanford University.
Bug finding is always appreciated, as well as corner cases where your protein
won't dabble.

- Issue tracker: https://github.com/drorlab/dabble/issues
- Source code: https://github.com/drorlab/dabble/project

License

-------
This projects is licensed under the GPLv2 license.


API Documentation
-----------------

.. toctree::
   :maxdepth: 2

   Dabble.rst
   DabbleParam.rst
   DabbleUtils.rst

* :ref:`search`

