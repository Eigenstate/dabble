## Dabble ##

*Dabble* is a Python program that facilitates the insertion of a protein system
into a pre-equilibrated membrane. It incorporates the atomselection language
from VMD, which is compiled into a Python module and imported into the program.

To use, just run dabble from this directory. It will pull in the relevant
dependencies and launch.

## Basic Usage ##

### The protein ###
Dabble takes as input a pdb or mae file containing a protein of interest.
It infers the file type from the file name.

      -i <file.mae/file.pdb>

The protein will be centered at the origin and the membrane will be added 
around it on the XY-plane. If your protein is not oriented right or you need
to move it, there are several ways to fix it.

First, you can provide an OPM structure to which the protein should be aligned.
Note that this protein must have the same number of atoms as the input protein.
You can specify which atoms should be aligned as well using VMD atom selection
language, which defaults to atoms that are "protein and backbone".

      --opm-pdb <opm.pdb> --opm-align "protein and backbone"

Alternatively, you can manually specify the orientation of the protein relative
to the membrane in terms of its angle to the membrane and z offset. The membrane
angle is the rotation of the membrane relative to the axis of the protein, in
degrees, as on the OPM website. Z offset is applied to the protein. The membrane
is always centered at the origin and on the XY plane.

      --move-solute <num> --membrane-rotation <num>


### The membrane ###

You can also specify a membrane to insert into. The membrane should be equilibrated,
and solvated with enough water on either side to fit your protein. It is best to
prepare your membrane with a lot of extra water, as it will be truncated to fit the
protein. The membrane should be oriented on the XY plane, and centered at the origin,
and in mae file format.

      -M <membrane.mae>

Dabble uses an equilibrated POPC membrane if you don't specify a membrane of your own.

### The output ###

Dabble is a one stop shop for all your parameterization needs, and supports output
in a variety of formats. 

However, if you are going to run a simulation, dabble can construct the membrane system
and then apply atom types and connectivity from your favourite molecular dynamics 
force field. Currently this only works if your favourite molecular dynamics force field
is CHARMM and you want to simulate in CHARMM or AMBER.

All intermediate files are saved, so if you want an AMBER prmtop, you will also get
the .mae, .pdb, and .psf files that were created along the way. These can be useful in
validating your structure and also save you the time of running dabble more than once.

#### Structure only ####
If you just want a structure, dabble can give you a mae, dms, or pdb file.

    -o <output.dms/output.mae/output.pdb>

The desired file format will be guessed from the extension, and no additional input
is needed.

#### CHARMM-ready ####

Dabble interfaces with psfgen to create a protein structure file (psf) and coordinate 
file (pdb) describing the system with CHARMM atom types. Dabble does all of the heavy
lifting, most of which is dealing with quirks of psfgen and handling translation from
pdb atom names/types to CHARMM ones. 

By default, the charmm36 atom names, types, and topologies aare used.
You will be prompted for additional or alternative rtf files where you can define
ligands or use a different version of the force field.

Dabble will attempt to match up atom names from the structure to the topologies, but
will occasionally fail. If this happens, you will be presented with a list of unmatched
atoms for each problem residue and a list of atoms in your file for that residue that could
not be matched. Usually it is easy to identify the correspondence: 1H -> H1, etc. However
if it is not obvious, please check the charmm rtf file and your structure to be sure that
the atoms are being done right.

To get a pdb and psf file,

    -o <output.psf>

And you will get two files, output.psf and output.pdb.

#### AMBER-ready with CHARMM parameters ####

No, you can't use AMBER parameters yet, sorry. But if you want to simulate in AMBER 
with CHARMM parameters, dabble can help you out. It will obtain a psf file using the same
method described above, and then will use the chamber function of the ParmEd API (by
Jason Swails) to produce AMBER input files.

Please be patient. The call to chamber often takes some time.

To get AMBER input files:

    -o <output.prmtop>

And you will get two files, output.prmtop and output.inpcrd, ready to run in AMBER.
_Please note that these files will not view correctly in VMD._ There will be a complaint
about the CTITLE record and no bonds will appear. This is due to VMD incorrectly parsing
the prmtop, not due to any errors in the process. I recommend loading the intermediate psf
file instead of the prmtop to check the final structure.


## More advanced usage, by example ##

## *NEW* Custom membranes, or no membranes ##

*"I want a system with just water, no membrane*

Lucky for you, there is new functionality to support this for TIP3 waters.

    -M TIP3

will use Dabble's pre-equilibrated box of TIP3 water as a solvent, with no
membrane. Note that this box is pretty small, so you will see some tiling effects.

*"I want to provide my own membrane*

There are a few steps to get a membrane in the correct format for dabbling.
The easiest way to get a membrane is to use [CHARMM-GUI.](http://www.charmm-gui.org/?doc=input/membrane)

1. Specify "membrane-only" system. Then build a membrane of desired composition. An
XY dimension of 30-50 will produce good results, although dabble can tile as necessary if
the provided membrane is too small. Continue to the next step in CHARMM-GUI.

2. No need to place ions as dabble will do this for you. Continue to the next step in CHARMM-GUI.

3. Stop following the "Assemble Components" step. Download the assembled psf and crd. NOT the pdb!
(step5_assembly).

4. Use the provided script to set a periodic box on the system. Dabble needs to know this informatoin
in order to correctly tile it. The script will prompt you for the directory in which the psf and
crd are saved, and will output a step5_assembly_dabble.mae file.
    
    convert_step5_to_dabble.py 

5. Rename your converted membrane so you don't forget what it is 
    mv step5_assembly_dabble.mae POPC_POPE_1-1.mae

5. Check the membrane looks and tiles correctly by visualizing it in VMD. If not, adjust the
periodic box as necessary.

6. You can now dabble it. Extra water will be added as necessary to accomodate your protein.
To make sure that your lipids are recognized correctly, you will probably have to provide an
atom selection for the correct residues.

    dabble.py -M POPC_POPE_1-1.mae --lipid-selection "resname POPC POPE" ...

## Salt concentrations ##

*"It added salt! I don't want salt!"*

You can specify the desired salt concentration, in M. The default is 0.150 M

    --salt-concentration 0.0

*"I need potassium, not sodium, for the cation."*

You can choose whether to use Na or K for the cation. The default is Na.

    --cation K

### Lipid and protein interactions protein ###

*"I have a protein with a palmitoylation, and I want the palmitoyl to be stuck in
the lipids"*

Use the "lipid friendly selection" option to specify atoms that are lipid friendly
and allowed to be much closer to the lipid atoms than protein atoms normally would be.
The selection is done using VMD's atom select syntax. Usually the "resname" option
works well.

    --lipid-friendly-sel "resname PCYS"

*"I have a cholesterol but lipid tails might get stuck in the rings"*

Use the lipid clash check option to specify atoms that are part of lipid residues with
rings that may run into other lipids. Cholesterol is really the prime example of this.

    --lipid-clash-check "resname CLOL"

*"I want the protein to be much farther from the lipids"*

You can specify the minimum distance between protein and lipid residues to move the
lipid either closer or farther from your protein. The default is 1.75 A.

    --lipid-dist 2.0


### Custom system dimensions ###

*"I need lots of water on the top and bottom of the protein, and dabble keeps
cutting it off too close"*

There are two ways of doing this. First, you can specify the desired absolute size of
your system, which will take priority over a buffer-based calculation. You can specify
only one dimension, in which case the other will be buffer based:

    --absolute-z 30.0

Or you can specify the entire size of the system:
    
    --absolute-dim 20.0,30.0

You could also change the buffer values, which can be easier if you don't know the 
exact size of your protein or want a guaranteed amount of water or lipid.

    --z-buffer-dist 10.0


*"I want a much larger area of membrane around the protein"*

You can specify the buffer distance around the protein that must consist of lipids.
The default is 20.0 A, which errs on the small side.

    --membrane-buffer-dist 30.0

Alternatively, you could specify the desired XY dimension of the final system:

    --absolute-xy 20.0

Or you can specify the entire size of the system:
    
    --absolute-dim 20.0,30.0


## Troubleshooting ##

*"I asked for a membrane system, but my protein ended up being just in water?'*

Short answer: Your initial protein is oriented wrong. Align it to the OPM 
structure and re-dabble.

Long answer: Dabble treats the Z coordinate of the protein as "truth", and will
only center it in the XY plane. Following insertion into the membrane, which
is centered at (0,0,0), extra atoms far from the protein are removed. If the Z
dimension of the protein is wrong, it can end up in the water away from the membrane,
which is then trimmed away.

*"I'm seeing waters with and without the TIP3 fake bond in my mae/pdb/dms. What's happening?"*

If you're simulating with TIP3 waters, this isn't a problem for simulation, as the
bond will be added in during parameterization steps (psf, prmtop generation, etc). 
What happened is your input membrane had waters with the bond explicitly defined, but
there was not enough water to solvate your protein with the desired buffer. More waters
are added, and dabble's reference water does not have the bond explicitly defined in
the mae topology. 


