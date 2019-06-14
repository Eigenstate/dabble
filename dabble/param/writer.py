"""
This module contains the MoleculeWriter class. It is used to apply
atom names from known topologies to the molecule by using a graph-based
representation of each molecule.

Author: Robin Betz

Copyright (C) 2019 Robin Betz
"""

# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330
# Boston, MA 02111-1307, USA.


from __future__ import print_function
import os
import logging

from abc import ABC, abstractmethod
from pkg_resources import resource_filename
from vmd import atomsel

logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class MoleculeWriter(ABC):
    """
    Represents a writer that will parameterize built systems into input files
    for some simulation program.

    Attributes:
        molid (int): VMD molecule ID to write
        output_prefix (str): Prefix for output file names. Appropriate suffix
            will be added depending on the file type.
    """

    #==========================================================================

    def __init__(self, molid, **kwargs):
        """
        Creates a molecule writer

        Args:
            molid (int): VMD molecule ID of system to write
            tmp_dir (str): Directory for temporary files. Defaults to "."
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            forcefield (str): Force field to use
            hmr (bool): If hydrogen masses should be repartitioned. Defaults
                to False.

            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default forcefield parameters
            debug_verbose (bool): Prints additional output, like from tleap.
        """
        self.molid = molid
        self.outprefix = ""
        self.matcher = None

        # Set default options
        self.tmp_dir = kwargs.get("tmp_dir", os.getcwd())
        self.lipid_sel = kwargs.get("lipid_sel", "lipid")
        self.debug = kwargs.get("debug_verbose", False)
        self.override = kwargs.get("override_defaults", False)
        self.hmr = kwargs.get("hmr", False)

        self.extra_topos = kwargs.get("extra_topos")
        self.extra_params = kwargs.get("extra_params")
        # Handle None from argparse in command line invocation
        if self.extra_topos is None: self.extra_topos = []
        if self.extra_params is None: self.extra_params = []

        # Handle None from argparse in command line invocation
        if self.extra_topos is None: self.extra_topos = []
        if self.extra_params is None: self.extra_params = []

    #==========================================================================

    def _apply_naming_dictionary(self, atomnames, resnames, verbose=False):
        """
        Applies the atom names from a matcher.

        Args:
            atomnames (dict int->str): Atom index to atom name
            resnames (dict int->str or str): Atom index to residue name,
                or just residue name for all atoms
            verbose (bool): If renamings should be printed out

        Raises:
            ValueError: If indices in atomnames and resnames dictionary
                differ
        """
        if isinstance(resnames, dict) and \
           set(resnames.keys()) != set(atomnames.keys()):
            raise ValueError("Invalid matching dictionary for resnames '%s'"
                             % resnames.keys())

        for idx, name in atomnames.items():
            atom = atomsel("index %s" % idx)
            # TODO: Name short circuiting here??
            if verbose and atom.name[0] != name:
                print("Renaming %s:%s: %s -> %s" % (atom.resname[0],
                                                    atom.resid[0], atom.name[0],
                                                    name))
            atom.name = name

            if isinstance(resnames, dict):
                atom.resname = resnames[idx]
            else:
                atom.resname = resnames

    #==========================================================================

    @abstractmethod
    def write(self, filename):
        pass

    @abstractmethod
    def get_topologies(forcefield):
        pass

    @abstractmethod
    def get_parameters(forcefield):
        pass

    #==========================================================================
    #                              Static methods
    #==========================================================================

    @staticmethod
    def get_pdb_line(atom, index, resindex, hetatom=False):
        """
        Get the PDB-formatted line corresponding to this atom, with a newline
        at the end.

        Args:
            atom (VMD atomsel): Atom selected
            index (int): Index in PDB file.
            resindex (int): Residue number in PDB file
            hetatom (bool): If this is part of a non-standard residue

        Returns:
            (str) PDB file line for this atom
        """
        if len(atom) != 1:
            raise ValueError("PDB entry selection must be only one atom")

        ins = atom.insertion[0]
        result = "%-6s%5d %-5s%-4s%c%4d%c  %8.3f%8.3f%8.3f%6.2f%6.2f" \
                 "     %-4s%2s\n" % ("HETATM" if hetatom else "ATOM",
                                     index,
                                     atom.name[0],
                                     atom.resname[0],
                                     atom.chain[0],
                                     resindex,
                                     ins if len(ins) else " ",
                                     atom.x[0],
                                     atom.y[0],
                                     atom.z[0],
                                     0.0, 0.0,
                                     atom.segname[0],
                                     atom.element[0])
        return result

    #==========================================================================

    @staticmethod
    def _get_forcefield_path(filename):
        """
        Gets the absolute path to a given bundled forcefield file
        """
        rn = resource_filename(__name__, os.path.join("parameters", filename))
        return os.path.abspath(rn)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
