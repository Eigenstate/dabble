# Sets appropriate environment variables
# and run
import sys, os
os.environ['DABBLEDIR'] = os.path.abspath(os.path.dirname(sys.argv[0]))+"/.."
os.environ['VMDDIR'] = os.environ['DABBLEDIR'] + "/vmd-1.9.2-python"
if os.environ['VMDDIR'] not in sys.path :
  sys.path.append( os.environ['VMDDIR'] )
  sys.path.append( os.environ['VMDDIR'] + "/scripts/python" )
  sys.path.append( os.path.abspath(os.path.dirname(sys.argv[0])) + "/src" )

import dabble
dabble.dabble(sys.argv)

