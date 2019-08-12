.. highlight:: bash
.. _command_line:

Command-line Options
====================

Dabble has many command line options. See :ref:`examples` for typical usage.

The command-line tool is ``dabble``. Run with no arguments
to see a list of all arguments, or for help::

   dabble -h

Input and output
----------------

Dabble takes as input a ``pdb``, ``mae``, or ``dms`` file containing a prepared protein
of interest::

   -i <protein file name>

"Prepared" in this context means hydrogens are added, caps (if desired) are added
to the end of protein chains, and missing atoms, residues, or loops have been filled in.
Dabble will not add or remove atoms from the input protein.

Specify the output file name with::

   -o <output>

The format of the output will be inferred by file extension, with the extension
referring to that of the topology file. Multiple output files will be produced
with different extensions, such as , or you can specify
it manually with the ``-format`` argument::

   --format [amber|charmm|desmond|gromacs|lammps|mae|pdb]

.. note::

   Sometimes, Dabble will produce multiple output format files, depending on
   the requested forcefield and format. For example, requesting a system be
   dabbled with the charmm forcefield with output in amber format will produce
   a ``.psf``, ``.pdb`` CHARMM-format files, along with ``.prmtop`` and
   ``.inpcrd`` AMBER-format files. The table only shows the files you are
   guaranteed to get regardless of chosen forcefield.

The relationship between file format, file extension, and format string is:

.. list-table::
   :header-rows: 1
   :widths: 10 5 40 5

   * - Format
     - Extension
     - Supported codes
     - Produced files
   * - amber
     - ``.prmtop``
     - AMBER, NAMD, OpenMM, ACEMD
     - ``.prmtop``, ``.inpcrd``
   * - charmm
     - ``.psf``
     - CHARMM, NAMD, OpenMM, ACEMD
     - ``.psf``, ``.pdb``
   * - desmond
     - ``.dms``
     - Anton (via conversion), Desmond
     - ``.dms``
   * - gromacs
     - ``.gro``
     - GROMACS, OpenMM
     - ``.top``, ``.gro``
   * - lammps
     - ``.dat``
     - LAMMPS
     - ``.dat``, ``.psf``, ``.pdb``
   * - mae
     - ``.mae``
     - N/A - unparameterized
     - ``.mae``
   * - pdb
     - ``.pdb``
     - N/A - unparameterized
     - ``.pdb``


.. note::

   If building a system, Dabble will always produce a ``.mae`` file with the
   build, unparameterized system. This can be useful if testing multiple
   forcefields.

Dabble will never overwrite any existing files with its output. If you wish
to disable this functionality, you can enable overwriting with:::

   -O

Or ``--overwrite``.

Builder options
---------------

Geometry specification
~~~~~~~~~~~~~~~~~~~~~~

Dabble builds the system by inserting the protein into a membrane, then adding
water and ions. The exact dimensions of the membrane, and amount of water,
can all be specified.

I recommend using a buffer-based specification, where water and lipid are placed
so that there is the specified amount of padding in the X, Y, and Z dimensions.

The buffer distance is from the edge of the protein to the edge of the periodic
box, so the distance between protein images is double this value.::

    --membrane-buffer-dist 17.5 --water-buffer 10.0

For water-only systems, no need to specify the membrane buffer::

    -M none --water-buffer 15.0

The X and Y dimensions of the protein are handled separately using the
buffer-based calculation, so the resulting system may not be square.  If
instead you want to manually specify your system size, you can use the
following options. This will produce a system where the periodic box is 50 x 50
x 100 A.::

    --absolute-x 50.0 --absolute-y 50.0 --absolute-z 100.0

You can mix the buffer-based and absolute system size specifications. Any
absolute dimensions will take priority.::

    --membrane-buffer-dist 20.0 --absolute-z 50.0

Input system orientation
~~~~~~~~~~~~~~~~~~~~~~~~

The input system will be centered in the XY plane and the membrane will be
added around it on the XY plane. Water will be added in the Z dimension.  If
your protein or input molecule(s) is not oriented right or you need to move it,
there are several ways to fix it.

The simplest is to pre-align to the `Orientation of Proteins in Membranes
<http://opm.phar.umich.edu/>`_ database when preparing the protein for
dabbling.

Alternatively, Dabble can attempt to perform an alignment. Note that due to the
functionality available in the vmd-python backend, the alignment must be done
to equal numbers of atoms in the reference and target proteins. You can specify
which atoms should be aligned using a VMD atom selection language, which
defaults to ``protein and backbone``.::

    --opm-pdb <opm pdb file> --opm-align "protein and backbone"

You can also manually specify the orientation of the protein relative to the
membrane in terms of its angle to the membrane and Z offset. The membrane angle
is the rotation of the membrane relative to the axis of the protein, in
degrees, as on the OPM website. Z offset is applied directly to the protein's
coordinates. The membrane is always centered in the XY plane.::

    --move-solute <z offset> --membrane-rotation <degrees>


Lipid membrane
~~~~~~~~~~~~~~

Dabble can add a lipid membrane to your input structure. If no membrane is
specified, Dabble will use a POPC membrane. However, you may also specify your
own. The membrane should be equilibrated, and can include any amount of water
in the +- Z direction, as Dabble will trim excess. If there is insufficient
water to solvate your protein, Dabble will add more, but it will require
equilibration.

To build your own membranes for Dabble using CHARMM-GUI, please read
:ref:`membranes`

.. todo::

   Support for direct loading of CHARMM-GUI membranes in PDB format will
   be added in the future.

The membrane should be oriented on the XY plane and in mae file format.::

    -M <membrane mae file>

If you don't want a membrane, Dabble solvate the system in just water. The
water model that will be used is specified elsewhere.::

    -M none

Dabble tries to delete atoms from the membrane that run into the input structure.
Its default logic works well, but you can provide your own selection strings
if you have an unusual membrane. The defaults are shown here.

Sometimes Dabble may not recognize your custom membrane as being composed of
lipids. If this is the case, you can manually specify an atom selection for the
lipid residues. The default value will pull out the following resnames: ``DLPE
DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE POPG POPS``. The following selection
should correspond to the lipids in the input membrane:::

   --lipid-selection "resname DOPC"

You can also specify the minimum distance between protein and lipid residues
to move the lipid either closer or farther from the protein. The default
value is 1.75 A.::

    --lipid-dist 2.0

Lipid or integral membrane molecules that have rings that may clash with lipids
should be represented by the lipid clash selection:::

   --lipid-clash-check "resname CLR CLOL"

Parts of the input system that are "lipid-friendly" and are allowed to be closer
to the protein are in the lipid friendly selection. This is useful if you have
a palmitoylation or other post-translational molecules on an input protein:::

   --lipid-friendly-sel "none"


Ions
~~~~

Dabble will add ions in the solvent to the desired salt concentration (defaults
to 0.150 M, which is approximately physiological). Then, anions or cations will
be deleted until the system is neutral.

Currently the supported cations are Na:sub:`+` and K:sub:`+`, with Na:sub:`+`
being the default.::

    --cation K

The default anion is Cl:sub:`-`::

    --anion Cl

To add ions just so the system is neutral:::

    --salt-concentration 0.0

.. todo::

   Cations and anions with charge more than than 1 have not been tested,
   and the system may not be neutral.


Parameterization options
------------------------

After building a system, Dabble can assign parameters from a force field and
produce parameterized files suitable for input to a simulation code.

Dabble interfaces with the appropriate parameterization program (`psfgen`,
`tleap`, or `pdb2gmx`). It is an expert user of the program, setting atom types
and handling covalent linkages. It can even detect and apply CHARMM patches!

Force field
~~~~~~~~~~~

Specify the desired force field:::

   --forcefield "charmm"

Currently supported values are:

.. list-table::
   :header-rows: 1
   :widths: 5 20

   * - Forcefield
     - Description
   * - ``amber``
     - Amber 14: ff14SB protein, GAFF2 small molecule, lipid14 lipids
   * - ``charmm``
     - CHARMM36m, July 2018 update
   * - ``opls``
     - OPLS AA/M, 2001 amino acid dihedrals

.. todo::

   Lipids are unsupported with the OPLS AA/M forcefield due to uncertainty
   about parameter files. If you are an expert OPLS user, please contact the
   developers to help.


Water model
~~~~~~~~~~~

You can also select which water model to use:::

   --water-model "tip3"

Currently supported values are:

.. list-table::
   :header-rows: 1
   :widths: 5 20

   * - Model
     - Description
   * - ``tip3``
     - TIP3 model, from W.L. Jorgensen, J.Chandrasekhar, J.D. Madura;
       R.W. Impey, M.L. Klein; Comparison of simple potential functions
       for simulating liquid water; J. Chem. Phys. 79 926-935 (1983)
   * - ``tip4e``
     - TIP4P-Ewald, from H.W. Horn, W.C Swope, J.W. Pitera, J.D. Madura,
       T.J. Dick, G.L. Hura, T. Head-Gordon; J. Chem. Phys.
       20: 9665-9678 (2004)
   * - ``spce``
     - SPC/E model, from H.J.C. Berendsen, J. R. Grigera,
       T. P. Straatsma; The Missing Term in Effective Pair
       Potentials; J. Phys. Chem 1987, 91, 6269-6271 (1987)

.. todo::

   TIP4P/Ew will not work with AMBER-format ``.prmtop`` output using the
   CHARMM forcefield, most likely due to a bug in ParmEd. This is under
   investigation.

Additional parameter files
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have ligands or other residues that are not represented in the default
force field, you must provide both a topology and parameter file for these
residues. See :ref:`tutorials` for how to obtain these files.

Pass topology files with the `-top` flag, and parameter files with the `-prm`
flag. You can provide these flags multiple times:::

   -ff amber -top ligand1.off -par ligand1.frcmod -top ligand2.off -parm ligand2.frcmod

For AMBER, Dabble can also parse `.leaprc` files that can load multiple other
topology and parameter files. These are considered a topology file but can
be also passed as a parameter file:::

   -ff amber -top all_ligands.leaprc

Some formats, like CHARMM, will combine the topology and parameter files into
one. You should pass this file as both:::

   -ff charmm -top ligand.str -par ligand.str


Custom forcefields
~~~~~~~~~~~~~~~~~~

If you want to use your own force field files and not have Dabble load any for
you, specify the forcefield flag as usual and pass ``--override-defaults``,
and all your parameter files. 

For example, to use only the old GAFF forcefield in AMBER:::

   -ff amber -par $AMBERHOME/dat/leap/cmd/leaprc.gaff --override-defaults


Hydrogen mass repartitioning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run your simulations in AMBER with timesteps up to 4 fs, Dabble can use the
ParmEd API to conduct hydrogen mass repartitioning. This only works when
requesting output in AMBER formats.::

    --format amber --hmr


Debug options
-------------

Dabble writes temporary files to a directory that it tries to clean up when
it's done. To retain this directory in a specific location:::

   --tmp-dir /tmp/dabble/

And to be more verbose::

   --verbose


