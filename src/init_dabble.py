# Sets appropriate environment variables
# and run
import sys, os
os.environ['DABBLEDIR'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/.."

# Put vmd.so loadable module in the path
os.environ['VMDDIR'] = os.environ['DABBLEDIR'] + "/vmd-1.9.2-python"
if os.environ['VMDDIR'] not in sys.path :
  sys.path.append( os.environ['VMDDIR'] )
  sys.path.append( os.environ['VMDDIR'] + "/scripts/python" )
  sys.path.append( os.path.abspath(os.path.dirname(sys.argv[0])) + "/src" )

# Put the parmed installation in the path
sys.path.append( os.environ['DABBLEDIR'] + "parmed/src/ParmedTools" )
sys.path.append( os.environ['DABBLEDIR'] + "parmed/src/chemistry" )
sys.path.append( os.environ['DABBLEDIR'] + "parmed/src/chemistry/amber" )

import dabble
dabble.dabble(sys.argv)

