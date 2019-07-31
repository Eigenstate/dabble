"""
This module contains the LammpsWriter class. It invokes the appropriate
parameterizer then converts the results to input format suitable for simulation
in LAMMPS.

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

from vmd import evaltcl, molecule

from dabble import DabbleError
from dabble.param import (AmberWriter, CharmmWriter, MoleculeWriter)

logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class LammpsWriter(MoleculeWriter):
    """
    Represents a writer that will parameterize built systems into input files
    for input to LAMMPS.

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
            forcefield (str): Force field to use
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            hmr (bool): If hydrogen masses should be repartitioned. Defaults
                to False.
            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default charmm parameters
            debug_verbose (bool): Prints additional output, like from tleap.
        """
        super(LammpsWriter, self).__init__(molid, **kwargs)

        self.forcefield = kwargs.get("forcefield", "charmm")
        self.hmr = kwargs.get("hmr", False)

        # Our get_ method handles all forcefield combinations
        self.topologies = self.get_topologies(self.forcefield)
        self.parameters = self.get_parameters(self.forcefield)

        # Handle override
        if self.override:
            self.topologies = []
            self.parameters = []

        # Now extra topologies and parameters
        self.topologies.extend(self.extra_topos)
        self.parameters.extend(self.extra_params)

        self.matcher = None

    #==========================================================================

    def write(self, filename):
        """
        Writes the parameter and topology files.

        Args:
            filename (str): File name to write. Gromacs suffix will be added.
        """
        self.outprefix = filename

        # Charmm forcefield
        if "charmm" in self.forcefield or "opls" in self.forcefield:
            psfgen = CharmmWriter(molid=self.molid,
                                  tmp_dir=self.tmp_dir,
                                  lipid_sel=self.lipid_sel,
                                  forcefield=self.forcefield,
                                  water_model=self.water_model,
                                  hmr=self.hmr,
                                  extra_topos=self.extra_topos,
                                  extra_params=self.extra_params,
                                  override_defaults=self.override)
            print("Writing intermediate psf")
            psfgen.write(self.outprefix)
            self._psf_to_lammps()

        elif "amber" in self.forcefield:
            prmgen = AmberWriter(molid=self.molid,
                                 tmp_dir=self.tmp_dir,
                                 forcefield=self.forcefield,
                                 water_model=self.water_model,
                                 hmr=self.hmr,
                                 lipid_sel=self.lipid_sel,
                                 extra_topos=self.extra_topos,
                                 extra_params=self.extra_params,
                                 override_defaults=self.override)
            print("Writing intermediate prmtop")
            prmgen.write(self.outprefix)
            self._prmtop_to_lammps()

        else:
            # Currently unsupported
            raise DabbleError("Forcefield '%s' not supported for lammps"
                              % self.forcefield)

    #==========================================================================

    def _prmtop_to_lammps(self):

        # Load prmtop, save a psf, then do _psf_to_lammps on that
        m = molecule.load("prmtop", self.outprefix + ".prmtop",
                          "rst7", self.outprefix + ".inpcrd")
        molecule.write(m, "psf", self.outprefix + ".psf")
        molecule.write(m, "pdb", self.outprefix + ".pdb")

        self._psf_to_lammps()


    #==========================================================================

    def _psf_to_lammps(self):

        # Write a temporary file for a topotools tcl script
        f, temp = tempfile.mkstemp(prefix="topo", suffix=".tcl",
                                    dir=self.tmp_dir, text=True)
        with os.fdopen(f, 'w') as fh:
            fh.write("package require topotools 1.6\n")
            fh.write("mol load psf %s.psf pdb %s.pdb\n"
                     % (self.outprefix, self.outprefix))
            fh.write("topo writelammpsdata %s.dat full" % self.outprefix)

        # Run the pbctools tcl script
        output = evaltcl("play %s" % temp)
        print(output)

    #==========================================================================
    #                              Static methods
    #==========================================================================

    @classmethod
    def get_topologies(cls, forcefield):
        """
        Gets topologies depending on the forcefield
        """

        # Amber, Charmm, and OPLS handled by conversion
        if forcefield == "charmm":
            return CharmmWriter.get_topologies(forcefield)

        if forcefield == "amber":
            return AmberWriter.get_topologies(forcefield)

        if forcefield == "opls":
            return CharmmWriter.get_topologies(forcefield)

        raise DabbleError("Unsupported forcefield %s" % forcefield)

    #==========================================================================

    @classmethod
    def get_parameters(cls, forcefield):
        """
        Gets the parameters depending on the forcefield
        """

        # Amber and Charmm handled by converstion
        if forcefield == "charmm":
            return CharmmWriter.get_parameters(forcefield)

        if forcefield == "amber":
            return AmberWriter.get_parameters(forcefield)

        if forcefield == "opls":
            return CharmmWriter.get_parameters(forcefield)

        raise DabbleError("Unsupported forcefield %s" % forcefield)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
