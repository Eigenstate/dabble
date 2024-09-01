.. highlight:: bash
.. _installation:

Installation
============

The preferred installation for Dabble is with ``conda``.
This will install all dependencies automatically::

   conda install -c rbetz dabble

Anaconda (or conda) is a package manager for Python that simplifies managing
multiple Python versions and package dependencies. If you are new to scientific
Python, read more about the `Anaconda Python distribution
<https://www.anaconda.com/distribution/>`_.


From Source
-----------

Dabble depends on the following Python packages:

- `vmd-python <https://vmd.robinbetz.com/>`_
- `psfgen <https://psfgen.robinbetz.com/>`_
- `networkx <https://networkx.github.io/>`_
- `ParmEd <https://parmed.github.io/ParmEd/html/index.html/>`_

You will also need an installation of AmberTools for AMBER-format
output, and Gromacs from Gromacs-format output.

With all dependencies installed, you can install Dabble from the
`github repository <https://github.com/Eigenstate/dabble>`_ as follows::


   git clone git@github.com:Eigenstate/dabble.git
   cd dabble/
   python setup.py install



