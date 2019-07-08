"""
Parameterizes molecules for molecular dynamics simulations
"""

__version__ = '2.7.10'
__author__ = 'Robin Betz'

# Currently supported forcefields and information
supported_forcefields = {
    "charmm": "CHARMM36m, July 2018 update",
    "amber": "AMBER 14",
    "opls": "OPLS AA/M, 2001 aminoacid dihedrals",
}

from dabble.param.moleculematcher import MoleculeMatcher
from dabble.param.ambermatcher import AmberMatcher
from dabble.param.charmmmatcher import CharmmMatcher, Patch
from dabble.param.gromacsmatcher import GromacsMatcher

from dabble.param.writer import MoleculeWriter
from dabble.param.amber import AmberWriter
from dabble.param.charmm import CharmmWriter
from dabble.param.gromacs import GromacsWriter
from dabble.param.lammps import LammpsWriter

