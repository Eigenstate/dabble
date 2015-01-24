## Dabble ##

*Dabble* is a Python program that facilitates the insertion of a protein system
into a pre-equilibrated membrane. It incorporates the atomselection language
from VMD, which is compiled into a Python module and imported into the program.

To use, just run dabble from this directory. It will pull in the relevant
dependencies and launch.

The following command-line flags are mandatory:
- -i/--input <solute filename>    Path to the file containing the system to solvate.
  Currently only -mae format is supported.
- -o/--output <output filename>    Desired output name. Currently supported formats are
mae and pdb. The format to use is inferred from the extension.
- -M/--membrane-system <membrane filename>  Path to the file containing a pre-equilibrated
  membrane in which to insert the solute. The membrane will be tiled in the
  x and y directions as needed to accomodate the solute. 

The following command-line flags are commonly helpful:
- -c/--cation [Na/K]    The cation to use
- -s/--salt-concentration <num> The molar concentration of ions to have
- -f/--lipid-friendly-sel "atomsel" Atom selection for parts of the protein that
  are allowed to clash with the lipids (such as lipid-modified residues, etc)
- -C/--lipid-clash-check "atomsel" Atom selection for lipids with rings that may
  run into other lipids (such as cholesterol)

For a description of all command line options, including options that allow more
precise specification of system size, please run.
    
      dabble -h

Also ask me (Robin B.) for help. 

