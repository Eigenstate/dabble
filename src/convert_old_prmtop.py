#!/share/PI/rondror/software/miniconda/bin/python

# Converts a CHARMM format prmtop to a vmd compatible one
import sys
from ParmedTools import parmout
from chemistry.amber import AmberParm

if len(sys.argv) is not 3 :
  print("USAGE: convert_old_prmtop.py <old prmtop> <new prmtop>")
  quit(1)

old_prmtop = sys.argv[1];
new_prmtop = sys.argv[2];

parm = AmberParm(old_prmtop)
write = parmout(parm, "%s /dev/null vmd" % new_prmtop)
write.execute()
