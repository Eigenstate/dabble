"""
This module contains the GromacsWriter class. It is used to apply
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
import tempfile

from dabble import DabbleError
from dabble.param import AmberWriter, CharmmWriter
from parmed.formats.registry import load_file
from parmed.gromacs import GromacsGroFile, GromacsTopologyFile
from vmd import atomsel, evaltcl, molecule

logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class GromacsWriter(object):
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
            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default charmm parameters
            debug_verbose (bool): Prints additional output, like from tleap.
        """
        self.molid = molid
        self.forcefield = kwargs.get("forcefield", "charmm")
        if self.forcefield not in ["amber", "charmm"]:
            raise DabbleError("Unsupported forcefield: %s" % self.forcefield)

        self.tmp_dir = kwargs.get("tmp_dir", ".")
        self.lipid_sel = kwargs.get("lipid_sel", "lipid")
        self.extra_topos = kwargs.get("extra_topos", [])
        self.extra_params = kwargs.get("extra_params", [])
        self.debug = kwargs.get("debug_verbose", False)
        self.override_defaults = kwargs.get("override_defaults", False)
        self.outprefix = ""

    #==========================================================================

    def write(self, filename):
        """
        Writes the parameter and topology files.

        Args:
            filename (str): File name to write. Gromacs suffix will be added.
        """
        self.outprefix = filename

        # Charmm forcefield
        if "charmm" in self.forcefield:
            # Topologies used will be found and returned by CharmmWriter
            self.topologies = CharmmWriter.get_topologies(self.forcefield)
            self.parameters = CharmmWriter.get_parameters(self.forcefield)

            psfgen = CharmmWriter(molid=self.molid,
                                  tmp_dir=self.tmp_dir,
                                  lipid_sel=self.lipid_sel,
                                  extra_topos=self.extra_topos,
                                  override_defaults=self.override_defaults)
            psfgen.write(self.outprefix)
            self._psf_to_gromacs()

        elif "amber" in self.forcefield:
            prmgen = AmberWriter(molid=self.molid,
                                 tmp_dir=self.tmp_dir,
                                 forcefield=self.forcefield,
                                 lipid_sel=self.lipid_sel,
                                 extra_topos=self.extra_topos,
                                 extra_params=self.extra_params)
            print("Writing intermediate prmtop")
            prmgen.write(self.outprefix)
            self._amber_to_gromacs()

    #==========================================================================

    def _amber_to_gromacs(self):

        # Load prmtop and inpcrd as a Structure
        parmstruct = load_file(self.outprefix + ".prmtop",
                               xyz=self.outprefix + ".inpcrd",
                               structure=True)

        # Save .gro coordinate file
        GromacsGroFile.write(struct=parmstruct,
                             dest=self.outprefix + ".gro")

        # Save .top topology and parameter file
        grotop = GromacsTopologyFile.from_structure(parmstruct, copy=False)
        grotop.write(dest=self.outprefix + ".top",
                     parameters="inline")

    #==========================================================================

    def _psf_to_gromacs(self):

        # Save the .gro coordinate file with VMD
        mid = molecule.load("psf", self.outprefix + ".psf",
                            "pdb", self.outprefix + ".pdb")
        atomsel("all", mid).write("gro", self.outprefix + ".gro")
        # TODO I think this doesn't write bonds!!!
        molecule.delete(mid)

        # Write a temporary file for a topotools tcl script
        f, temp = tempfile.mkstemp(prefix="topo", suffix=".tcl",
                                    dir=self.tmp_dir, text=True)
        with os.fdopen(f, 'w') as fh:
            fh.write("package require topotools 1.6\n")
            fh.write("mol load psf %s.psf pdb %s.pdb\n"
                     % (self.outprefix, self.outprefix))
            fh.write("topo writegmxtop %s.top [list %s ]\n"
                     % (self.outprefix, " ".join(self.parameters)))

        # Run the pbctools tcl script
        output = evaltcl("play %s" % temp)
        print(output)

    #==========================================================================

    # Satisfy inheritance requirement for these methods
    # I don't natively handle gromacs parameters so this doesn't do much
    def get_topologies(self, forcefield):
        return self.topologies

    def get_parameters(self, forcefield):
        return self.parameters

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
