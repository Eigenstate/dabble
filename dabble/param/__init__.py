"""
Parameterizes molecules for molecular dynamics simulations
"""

__version__ = '2.7.12'
__author__ = 'Robin Betz'

# Currently supported forcefields and information
supported_forcefields = {
    "charmm": "CHARMM36m, July 2018 update",
    "amber": "AMBER 14",
    "opls": "OPLS AA/M, 2001 aminoacid dihedrals",
}

supported_water_models = {
    "tip3": "TIP3 model, from W.L. Jorgensen, J.Chandrasekhar, J.D. Madura; "
            "R.W. Impey, M.L. Klein; Comparison of simple potential functions "
            "for simulating liquid water; J. Chem. Phys. 79 926-935 (1983).",
    "tip4p": "TIP4P-Ewald, from H.W. Horn, W.C Swope, J.W. Pitera, J.D. Madura,"
             " T.J. Dick, G.L. Hura, T. Head-Gordon; J. Chem. Phys. "
             "120: 9665-9678 (2004)",
    "spce": "SPC/E model, from H.J.C. Berendsen, J. R. Grigera, "
            "T. P. Straatsma; The Missing Term in Effective Pair "
            "Potentials; J. Phys. Chem 1987, 91, 6269-6271 (1987)",
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

