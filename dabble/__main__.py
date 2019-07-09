"""
Dabble, a membrane protein system builder

Author: Robin Betz

Copyright (C) 2019 Robin Betz

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330
Boston, MA 02111-1307, USA.

"""

from __future__ import print_function
import argparse
import os
import shutil
import signal
import sys
import tempfile
from dabble import VmdSilencer, DabbleBuilder, supported_formats
from dabble.param import supported_forcefields
from pkg_resources import resource_filename

__version__ = '2.7.11'
__author__ = 'Robin Betz'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Handle interrupts
def signal_handler(*args, **kwargs): # pylint: disable=unused-argument
    """ Catch signals """
    sys.stdout.write('\nInterrupted\n')
    sys.exit(1)

#==============================================================================

class DabbleTempDir(object): # pylint: disable=too-few-public-methods
    """
    Creates a destroyable temporary directory, but also allows it not
    to be destroyed, or it to be in a custom location.
    """
    def __init__(self, retain=False, path=None):
        if path is not None:
            if not os.path.isdir(path):
                os.mkdir(path)
            self.dir = path
            self.retain = True
        else:
            self.dir = tempfile.mkdtemp(prefix="dabble", dir=os.getcwd())
            self.retain = retain

    def __enter__(self):
        return self.dir

    def __exit__(self, types, value, traceback):
        if not self.retain:
            shutil.rmtree(self.dir, ignore_errors=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main(argv=None):

    WELCOME_SCREEN = '''
     ===============================================
    |               _      _      _                 |
    |             >(.)__ <(.)__ =(.)__              |
    |              (___/  (___/  (___/              |
    |                                               |
    |                    DABBLE            ______   |
    |                 _      _      _     /       \\ |
    |              __(.)< __(.)> __(.)=  <  beta! | |
    |              \\___)  \\___)  \\___)    \\_______/ |
    |                                               |
    |                 Robin Betz, 2019              |
    |               Stanford University             |
    | %s |
     ===============================================
    ''' % ('{0:^45}'.format("Version " + __version__))

    # pylint: disable=invalid-name
    parser = argparse.ArgumentParser(prog='dabble')

    group = parser.add_argument_group('Input and Output Files')
    group.add_argument('-i', '--input', dest='solute_filename',
                       metavar='<input>', type=str,
                       required=True,
                       help="Path to input protein or ligand file"
                      )
    group.add_argument('-o', '--output', dest='output_filename',
                       metavar='<output>', type=str,
                       required=True,
                       help="Name of output file. If -format argument is "
                       "supplied, the appropriate extension will be added. If "
                       "not, format will be inferred by the extension here. "
                       "Currently supported extensions: .pdb, .mae, .psf, .dms,"
                       " .prmtop"
                      )
    group.add_argument('-format', '--format', dest='format',
                       metavar='<format>', type=str,
                       default=None,
                       choices=supported_formats.keys(),
                       help="Format of output file. Supported: %s"
                       % ", ".join(supported_formats.keys())
                      )
    group.add_argument('-M', '--membrane-system', dest='membrane_system',
                       type=str, metavar='<solvent>',
                       default=resource_filename(__name__, "lipid_membranes/popc.mae"),
                       help="Path to pre-built membrane + solvent block. Must "
                       "be a .mae file. See documentation to create your own. "
                       "Defaults to a POPC membrane",
                       )
    group.add_argument('-O', '--overwrite', dest='overwrite',
                       action='store_true',
                       help="Overwrite existing files with requested output",
                      )

    group = parser.add_argument_group('Parameterization Options')
    group.add_argument('-ff', '--forcefield', dest='forcefield',
                       type=str, metavar='<forcefield>',
                       default="charmm",
                       choices=supported_forcefields.keys(),
                       required=True, action="store",
                       help="Force field to use for parameterization. Currently"
                            " supported values: %s"
                             % ", ".join(supported_forcefields.keys())
                      )
    group.add_argument('--hmr', dest='hmassrepartition',
                       default=False,
                       action='store_true',
                       help="Repartition Hydrogen masses to allow up to 4fs "
                       "time steps. Currently supported for AMBER (.prmtop) "
                       "output format only"
                      )
    group.add_argument('-top', '--topology', dest='extra_topos',
                       type=str, metavar='<topology file>',
                       default=None,
                       action='append',
                       help="Additional topology (rtf, off, lib, leaprc) file "
                       "to include in parameterization"
                      )
    group.add_argument('-par', '--parameters', dest='extra_params',
                       type=str, metavar='<parameter file>',
                       default=None,
                       action='append',
                       help="Additional parameter (prm, lib, frcmod) file to "
                       "include in parameterization"
                      )

    group = parser.add_argument_group('Lipid Membrane Options')
    group.add_argument('-L', '--lipid-selection', dest='lipid_sel',
                       type=str,
                       default='lipid or resname POPS POPG',
                       help="Atom selection string (VMD syntax) for the lipids "
                       "in the membrane. Defaults to 'lipid or resname POPS'"
                      )
    group.add_argument('-C', '--lipid-clash-check', dest='clash_lipids',
                       type=str,
                       default='resname CLR CLOL',
                       help="Atom selection string (VMD syntax) for lipids or "
                       "other integral membrane molecules with rings (i.e. "
                       "cholesterol) that might clash with other lipids. "
                       "Defaults to 'resname CLR CLOL'",
                      )
    group.add_argument('-f', '--lipid-friendly-sel', dest='lipid_friendly_sel',
                       type=str,
                       help="Atom selection string (VMD syntax) for parts of "
                       "the protein that are 'lipid-friendly' and should not be"
                       "considered when calculating which lipids are clashing "
                       "with the protein (i.e.: lipid tails, palmitoylations). "
                       "Defaults to no selection."
                      )

    group = parser.add_argument_group('Ion Options')
    group.add_argument('--cation',
                       default='Na', type=str,
                       help='Specify element of cation. Defaults to "Na"'
                      )
    group.add_argument('--anion',
                       default='Cl', type=str,
                       help='Specify element of anion. Defaults to "Cl"'
                      )
    group.add_argument('-s', '--salt-concentration', dest='salt_conc',
                       type=float,
                       default=0.150,
                       help="Salt concentration in final system, in M. Defaults"
                       " to physiological 0.150M NaCl"
                      )

    group = parser.add_argument_group('System Size Options')
    z_buffer_opts = group.add_mutually_exclusive_group()
    z_buffer_opts.add_argument('-w', '--water-buffer', dest='wat_buffer',
                               type=float,
                               default=20.0,
                               help="Buffer, in A, from each side of the "
                               "protein/solute to the edge of the periodic box."
                               " Defaults to a conservative 20.0 A"
                              )
    group.add_argument('-m', '--membrane-buffer-dist', dest='xy_buf', default=17.5,
                       type=float, help='membrane buffer distance from the protein to the '
                                        'box edge in the XY plane.'
                                        '[default: 17.5 angstroms]')
    group.add_argument('-d', '--lipid-dist', dest='lipid_dist',
                       default=1.75, type=float,
                       help='minimum distance from solute to lipid acyl group'
                       '[default: 1.75]')
    group.add_argument('--absolute-x', type=float, default=None,
                       dest='user_x', help='Specifies the x dimension. Takes '
                       'precedence over buffer-based calculation.')
    group.add_argument('--absolute-y', type=float, default=None,
                       dest='user_y', help='Specifies the y dimension. Takes '
                       'precedence over buffer-based calculation.')
    group.add_argument('--absolute-z', type=float, default=None,
                       dest='user_z', help='Specifies the z dimension. Takes '
                       'precedence over buffer-based calculation.')

    group = parser.add_argument_group('Orientation Options',
                                      'These options control how the input solute '
                                      'is oriented before inserting it into the '
                                      'solvent. Although it is recommended you '
                                      'pre-align the solute, these options are '
                                      'here for your convenience.')
    group.add_argument('--opm-pdb', dest='opm_pdb',
                       default=None, type=str,
                       help='oriented pdb file from OPM to align protein to'
                            '[default: None]')
    group.add_argument('--opm-align', dest='opm_align',
                       default='protein and backbone', type=str,
                       help='atomsel for OPM backbone atoms to align to'
                            '[default: protein and backbone]')
    group.add_argument('--move-solute', dest='z_move',
                       default=None, type=float,
                       help='value added to solute z coordinates'
                            '[default: 0]')
    group.add_argument('--membrane-rotation', dest='z_rotation',
                       default=None, type=float,
                       help='Membrane rotation relative to Z axis of protein, in '
                            'degrees. Use the number from OPM if you have it. '
                            '[default: 0]')

    group = parser.add_argument_group('Debug and Testing Options')
    group.add_argument('--tmp-dir', dest='tmp_dir', default=None)
    group.add_argument('--verbose', dest='debug_verbose', default=False,
                       action='store_true')

    print(WELCOME_SCREEN)
    print("\nCommand was:\n  %s\n" % " ".join([i for i in sys.argv]))
    #opts = parser.parse_args(sys.argv[1:])
    opts = parser.parse_args(argv)

# Make the temporary directory. Needs to be done now so there is somewhere
# to save the vmd output
    with DabbleTempDir(path=opts.tmp_dir, retain=opts.debug_verbose) as tdir:

        opts.tmp_dir = tdir # Needs to be defined in opts for builder to work
        soutput = sys.stdout if opts.debug_verbose else os.path.join(tdir,
                                                                     "vmd_output.txt")

        with VmdSilencer(output=soutput):
            signal.signal(signal.SIGINT, signal_handler)
            builder = DabbleBuilder(**vars(opts))
            builder.write()
            sys.stdout.flush()
            print("\nSuccess!")

if __name__ == "__main__":
    main()

