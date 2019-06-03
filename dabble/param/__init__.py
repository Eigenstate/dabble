"""
Parameterizes molecules for molecular dynamics simulations
"""

__version__ = '2.7.9'
__author__ = 'Robin Betz'

# Currently supported forcefields and information
supported_forcefields = {
    "charmm": "CHARMM36m, July 2018 update",
    "amber": "AMBER 14",
    "opls": "OPLS AA/M"
}

from dabble.param.moleculematcher import MoleculeMatcher
from dabble.param.charmmmatcher import CharmmMatcher, Patch
from dabble.param.ambermatcher import AmberMatcher
from dabble.param.writer import MoleculeWriter
from dabble.param.charmm import CharmmWriter
from dabble.param.amber import AmberWriter
from dabble.param.gromacs import GromacsWriter

