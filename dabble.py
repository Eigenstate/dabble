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
import signal
import sys
from Dabble import __version__, DabbleBuilder

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
|                 Robin Betz, 2015              |
|               Stanford University             |
| %s |
 ===============================================
''' % ('{0:^45}'.format("Version " + __version__))

# Handle interrupts
def signal_handler(*args, **kwargs): # pylint: disable=unused-argument
    """ Catch signals """
    sys.stdout.write('\nInterrupted\n')
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

# Log messages
def _make_logger(out, quiet=False):
    """
    Creates a logger that auto-flushes
    """
    def logger(msg): # pylint: disable=missing-docstring
        if not quiet:
            out.write(msg)
            out.flush()
    return logger

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
                   'prmtop (amber with charmm params)')
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
group.add_argument('--hmr', dest='hmassrepartition', default=False,
                   action='store_true', help='Repartition Hydrogen masses'
                   'to allow up to 4fs time steps. Currently amber only')
group.add_argument('-top', '--topology', default=None,
                    type=str, metavar='<topologies>', dest='extra_topos',
                    help='Additional topology (rtf or str) files to'
                    'include in parameterization, separated by commas')
group.add_argument('-par', '--parameters', default=None,
                   type=str, metavar='<parameters>', dest='extra_params',
                   help='Additional parameter (prm or str) files to'
                   'include in parameterization, separated by commas')

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
group.add_argument('-m', '--membrane-buffer-dist', dest='xy_buf', default=35.0,
                   type=float, help='buffer distance through the membrane.'
                   '[default: 35.0 angstroms]')
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
log = _make_logger(sys.stdout, opts.quiet)
log('\n\n')

builder = DabbleBuilder(**vars(opts)) # pylint: disable=star-args
builder.write()

