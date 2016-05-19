#!/usr/bin/env python
"""
Dabble, a membrane protein system builder

Author: Robin Betz

Copyright (C) 2015 Robin Betz

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
import os, sys
import signal
import tempfile

__version__ = '2.0.3'
__author__ = 'Robin Betz'
    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class VmdSilencer:
    """
    Toggles whether or not C extensions can write to stdout. Since
    VMD is the only C extension that does this, silences this extra
    output during the dabbling process. This is done in dabble.py not
    in the Dabble API because the user probably wants more info / it 
    can't hurt them if they're smart enough to use the API.

    Most of this is from:
    http://code.activestate.com/recipes/577564-context-manager-for-low-level-redirection-of-stdou/
    which is licensed under the MIT license.

    Attributes:
        output (file): Where to put the VMD ouptut
    """

    def __init__(self, output=os.devnull):
        self.outfile = output
        self.mode = 'w'
        
    #==========================================================================

    def __enter__(self):
        self.sys = sys
        # save previous stdout/stderr
        self.saved_stream = sys.__stdout__
        self.fd = self.saved_stream.fileno()
        self.saved_fd = os.dup(sys.stdout.fileno())
        sys.stdout.flush() # flush any pending output 

        # open surrogate files
        null_fd = open(self.outfile, self.mode)
        os.dup2(null_fd.fileno(), self.fd)

        self.null_stream = open(self.outfile, self.mode, 0)
        self.null_fd = self.null_stream.fileno()

        # overwrite file objects and low-level file descriptors
        os.dup2(self.null_fd, self.fd)

        sys.stdout = os.fdopen(self.saved_fd, 'w')
        sys.stdout.flush()

    #==========================================================================

    def __exit__(self, *args):
        sys = self.sys
        # flush any pending output
        sys.__stdout__.flush()
        # restore original streams and file descriptors
        os.dup2(self.saved_fd, self.fd)
        sys.stdout = self.saved_stream
        # clean up
        self.null_stream.close()
        os.close(self.saved_fd)
        return False

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Handle interrupts
def signal_handler(*args, **kwargs): # pylint: disable=unused-argument
    """ Catch signals """
    sys.stdout.write('\nInterrupted\n')
    sys.exit(1)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
|                 Robin Betz, 2016              |
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
                   help='Path to input protein or ligand file to build into '
                   'the system')
group.add_argument('-o', '--output', dest='output_filename',
                   metavar='<output>', type=str,
                   required=True,
                   help='Name of output file, format will be inferred by '
                   'extension. Currently supported: pdb, mae, psf (charmm), '
                   'prmtop (amber or charmm)')
group.add_argument('-M', '--membrane-system', dest='membrane_system',
                   type=str, metavar='<solvent>',
                   default="DEFAULT",
                   help='Path to solvent system (must be a mae file). Defaults '
                   'to a POPC membrane')
group.add_argument('-O', '--overwrite', dest='overwrite', action='store_true',
                   help='Overwrite output files, if found')
group.add_argument('-q', '--quiet', dest='quiet',
                   action='store_true', default=False)

group = parser.add_argument_group('Parameterization Options')
group.add_argument('-ff', '--forcefield', dest='forcefield',
                   type=str, metavar='<forcefield>', default="charmm",
                   choices=['amber','charmm'], action='store',
                   help="Force field to use for parameterization. Currently "
                        "supported: amber/charmm")
group.add_argument('--hmr', dest='hmassrepartition', default=False,
                   action='store_true', help='Repartition Hydrogen masses'
                   'to allow up to 4fs time steps. Currently prmtop output only')
group.add_argument('-top', '--topology', default=None, action='append',
                    type=str, metavar='<topologies>', dest='extra_topos',
                    help='Additional topology (rtf, off, lib) file to '
                    'include in parameterization')
group.add_argument('-par', '--parameters', default=None, action='append',
                   type=str, metavar='<parameters>', dest='extra_params',
                   help='Additional parameter (prm, lib, frcmod) file to '
                   'include in parameterization')
group.add_argument('-str', '--stream', default=None, action='append',
                   type=str, metavar='<streams>', dest='extra_streams',
                   help='Additional stream (str, leaprc) file to include in '
                   'parameterization')

group = parser.add_argument_group('Lipid Membrane Options')
group.add_argument('-L', '--lipid-selection', dest='lipid_sel',
                   default='lipid or resname POPS POPG', type=str,
                   help='atomsel for the lipids in the membrane [default: '
                   '"lipid or resname POPS"]')
group.add_argument('-C', '--lipid-clash-check', dest='clash_lipids',
                   help='Atomsel for lipids with rings (i.e. cholesterol) '
                   'that might clash with other lipids.')
group.add_argument('-f', '--lipid-friendly-sel', type=str,
                   dest='lipid_friendly_sel',
                   help='atomsel for parts of the protein that are '
                   '"lipid-friendly" and should not be used when calculating '
                   'which lipids are clashing with the protein (i.e.: lipid '
                   'tails, sidechains of peripheral membrane proteins)')
#
group = parser.add_argument_group('Ion Options')
group.add_argument('-c', '--cation',
                   default='Na', type=str,
                   help='specify cation "Na" or "K"'
                        '[default: "Na"]')
group.add_argument('-s', '--salt-concentration', dest='salt_conc',
                   default=0.150, type=float,
                   help='desired salt concentration.'
                        '[default: 0.150 M]')

group = parser.add_argument_group('System Size Options')
z_buffer_opts = group.add_mutually_exclusive_group()
z_buffer_opts.add_argument('-w', '--water-buffer', dest='wat_buffer', default=20.0,
                           type=float, help='water padding from each side of '
                           'protein [default 20.0 angstroms]')
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

print(WELCOME_SCREEN)
print("\nCommand was:\n  %s\n" % " ".join([i for i in sys.argv]))
opts = parser.parse_args(sys.argv[1:])

# Make the temporary directory. Needs to be done now so there is somewhere
# to save the vmd output
if not opts.tmp_dir:
    opts.tmp_dir = tempfile.mkdtemp(prefix='dabble', dir=os.getcwd())

with VmdSilencer(output=os.path.join(opts.tmp_dir,"vmd_output.txt")):

    signal.signal(signal.SIGINT, signal_handler)
    from Dabble import DabbleBuilder
    builder = DabbleBuilder(**vars(opts)) # pylint: disable=star-args
    builder.write()
    print("\nSuccess!")

