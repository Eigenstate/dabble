.. highlight:: bash
.. _command_line:

Command-line Usage
==================

The command-line tool is ``dabble``. Run with no arguments
to see a list of all arguments, or for help::

    dabble -h

The protein
-----------

Dabble takes as input a ``pdb``, ``mae``, or ``dms`` file containing a prepared protein
of interest:: 

    -i <protein file name>

"Prepared" in this context means hydrogens are added, caps (if desired) are added
to the end of protein chains, and missing atoms, residues, or loops have been filled in.
Dabble will not add or remove atoms from the input protein.

The protein will be centered in the XY plane and the membrane will be added around
it on the XY plane. Water will be added in the Z dimension.  
If your protein is not oriented right or you need to move it, there are several ways
to fix it.

The simplest is to pre-align to the `Orientation of Proteins in Membranes <http://opm.phar.umich.edu/>`_ database when preparing the protein for dabbling.

Alternatively, Dabble can attempt to perform an alignment. Note that due to the functionality available in the vmd-python backend, the alignment must be done to
equal numbers of atoms in the reference and target proteins. You can specify which
atoms should be aligned using a VMD atom selection language, which defaults to ``protein
and backbone``.::
    
    --opm-pdb <opm pdb file> --opm-align "protein and backbone"

You can also manually specify the orientation of the protein relative to the membrane
in terms of its angle to the membrane and Z offset. The membrane angle is the rotation
of the membrane relative to the axis of the protein, in degrees, as on the OPM website. Z offset is applied directly to the protein's coordinates. The membrane is always centered in the XY plane.::

    --move-solute <z offset> --membrane-rotation <degrees>

The membrane or solvent
-----------------------

If no membrane is specified, Dabble will use a POPC membrane. However, you may also specify your own. The membrane should be equilibrated, and can include any amount of water in the +- Z direction, as Dabble will trim excess. If there is insufficient water to solvate your protein, Dabble will add more, but it will require equilibration.

To build your own membranes for Dabble using CHARMM-GUI, please read 
:ref:`membranes`

The membrane should be oriented on the XY plane and in mae file format.::

    -M <membrane mae file>

If you don't want a membrane, Dabble can also solvate in TIP3 water with the
following option:::

    -M TIP3

Ions
----

Dabble will add ions in the solvent to the desired salt concentration (defaults
to 0.150 M, which is approximately physiological). Then, anions or cations will
be deleted until the system is neutral.

Currently the supported cations are Na:sub:`+` and K:sub:`+`, with Na:sub:`+`
being the default.::

    --cation K

To add ions just so the system is neutral:::

    --salt-concentration 0.0

Geometry specification
----------------------

Dabble builds the system by inserting the protein into a membrane, then adding
water and ions. The exact dimensions of the membrane, and amount of water,
can all be specified.

I recommend using a buffer-based specification, where water and lipid are placed
so that there is the specified amount of padding in the X, Y, and Z dimensions.

The buffer distance is from the edge of the protein to the edge of the periodic
box, so the distance between protein images is double this value.::

    --membrane-buffer-dist 17.5 --water-buffer 10.0

For water-only systems, no need to specify the membrane buffer::

    -M TIP3 --water-buffer 15.0

The X and Y dimensions of the protein are handled separately using the buffer-based
calculation, so the resulting system may not be square.
If instead you want to manually specify your system size, you can use the following
options. This will produce a system where the periodic box is 50 x 50 x 100 A.::

    --absolute-x 50.0 --absolute-y 50.0 --absolute-z 100.0

You can mix the buffer-based and absolute system size specifications. Any
absolute dimensions will take priority.::

    --membrane-buffer-dist 20.0 --absolute-z 50.0

Desired output
--------------

Dabble aims to be a one stop shop for all your system building and parameterization
needs, and supports output in a variety of formats, including common structure
files as well as parameterized input files for CHARMM, Desmond, and AMBER using
the CHARMM or AMBER force fields.

All intermediate files are saved, so if you want an AMBER prmtop, you will also get
the ``mae`` file made along the way. This can be useful in validating or 
visualizing your structure, and can save you the time of having to run Dabble
more than once.

Dabble will refuse to create output files if they would overwrite other files
with the same name. To allow overwriting to happen:::

    --overwrite

Structure only
~~~~~~~~~~~~~~

If you just want a built structure, Dabble can give you a ``mae``, ``dms``, or ``pdb``
file.::

    -o <output.dms/output.mae/output.pdb>

The file format will be guessed from the extension, and no additional input is needed.

Using CHARMM parameters
~~~~~~~~~~~~~~~~~~~~~~~

Dabble interfaces with `psfgen <http://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/>`_
to create a protein structure file (``psf``), and coordinate file (``pdb``) describing the system with CHARMM atom types and parameters.

Dabble does all the heavy lifting, most of which is dealing with the quirks of psfgen
and handling translation from pdb atom names and types to CHARMM ones. This matching functionality can even detect and apply patches! For more on how
this is done, please see :ref:`parameterization`.

By default, the charmm36 atom names, types, and topologies are used when you specify the ``-ff charmm`` flag indicating you want CHARMM parameters. If you wish to override this, pass the ``--no-default-ff`` flag and your own forcefield files.

Provide additional ``str``, ``rtf``/``top``, or ``par``/``prm`` files with the appropriate flag. Each flag may used multiple times in case you need to add multiple additional parameter sets.::

    -o output.psf -ff charmm -top ligand1.rtf -top ligand2.rtf -par ligands.prm -str ligand3.str

If you want to simulate in AMBER with CHARMM parameters, Dabble will invoke the `ParmEd
API <http://parmed.github.io/ParmEd/html/index.html>`_ to produce AMBER input files
with CHARMM parameters. Just request an AMBER format output topology file ``(prmtop)``
and specify CHARMM parameterization ``-ff charmm``. A coordinate file ``(inpcrd)``
will also be produced in AMBER format.::

    -o output.prmtop -ff charmm -str ligand.str

**NOTE:** These files may not view correctly in older versions of VMD. There will
be a complaint about the CTITLE record and no bond will appear. This is due to
VMD incorrectly parsing the prmtop and not any errors in the process. I have submitted
a patch to the VMD developers, but in the meantime I recommend loading the
intermediate ``psf`` file to check or visualize simulations.

AMBER parameters
~~~~~~~~~~~~~~~~

If you want to simulate in AMBER using AMBER parameters, Dabble interfaces with your local installation of `AmberTools <http://ambermd.org/#AmberTools>`_ to generate 
a topology file (``prmtop``) and coordinates (``inpcrd``) suitable for simulation.
You will need a local installation of AmberTools in the location specified by the
environment variable ``$AMBERHOME``.

Request AMBER parameterization with the ``-ff amber`` flag. By default,
the `ff14SB <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255>`_ protein parameter
set, `lipid14 <http://pubs.acs.org/doi/abs/10.1021/ct4010307>`_ lipid parameteres, 
and TIP3P water model, and GAFF2 small molecule parameters will be used. To override
these defaults, pass the ``--no-default-ff`` flag and your own forcefield leaprc files.

Provide parameter and residue definition files for ligands (``off`` or ``lib`` files)
using the ``-top`` flag, forcefield definition ``leaprc`` files with the ``-str`` flag,
and additional parameters ``frcmod`` with the ``-par`` flag. Dabble will look for
leaprc files in the ``$AMBERHOME/dat/leap/cmd`` directory as well as the current folder if no explicit path is provided.::

    -o output.prmtop -ff amber -top ligand.off -par ligand.frcmod -str leaprc.gaff

Less common options
-------------------

The following options are helpful, but not required to Dabble your system.

Hydrogen mass repartitioning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run your simulations in AMBER with timesteps up to 4 fs, Dabble can use the ParmEd API to conduct hydrogen mass repartitioning. This only works when requesting output
in AMBER formats.::
    
    -o output.prmtop --hmr

Lipid-protein interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Dabble inserts the protein into the membrane by tiling the membrane patch until
it is of appropriate size, overlapping the protein with it, and deleting any
lipids that are in the way. There are several ways to alter this behavior.

Use a "lipid friendly selection" to specify protein or ligand atoms that are
allowed to be much closer to the lipid atoms than protein atoms normally would
be. The selection is the VMD's atom select syntax. Usually the "resname" keyword
works well.::

    --lipid-friendly-sel "resname PCYS"

To specify atoms that are part of greasy or lipid residues with rings that
may run into other lipids, use the lipid clash check option. Cholesterol is
really the prime example of this.::

    --lipid-clash-check "resname CLOL"

You can also specify the minimum distance between protein and lipid residues
to move the lipid either closer or farther from the protein. The default
value is 1.75 A.::

    --lipid-dist 2.0

Sometimes Dabble may not recognize your custom membrane as being composed of
lipids. If this is the case, you can manually specify an atom selection for
the lipid residues. The default value will pull out the following resnames:
``DLPE DMPC DPPC GPC LPPC PALM PC PGCL POPC POPE POPG POPS``. Note that these
are resnames in the original input file provided to Dabble, as lipid selection
matters during the build phase.::

    --lipid-selection "resname DOPC"
