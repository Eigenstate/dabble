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
import subprocess

from parmed.formats.registry import load_file
from parmed.gromacs import GromacsGroFile, GromacsTopologyFile
from vmd import atomsel, evaltcl, molecule

from dabble import DabbleError
from dabble.param import (AmberWriter, CharmmWriter, MoleculeWriter,
                          GromacsMatcher)

logger = logging.getLogger(__name__) # pylint: disable=invalid-name

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class GromacsWriter(MoleculeWriter):
    """
    Represents a writer that will parameterize built systems into input files
    for input to GROMACS.

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
        super(GromacsWriter, self).__init__(molid, **kwargs)

        self.forcefield = kwargs.get("forcefield", "opls")
        self.hmr = kwargs.get("hmr", False)

        # Our get_ method handles all forcefield combinations
        self.topologies = self.get_topologies(self.forcefield)
        self.parameters = self.get_parameters(self.forcefield)
        self.ffdir = self.topologies[0]

        # Handle override
        if self.override:
            self.topologies = []
            self.parameters = []

        # Now extra topologies and parameters
        self.topologies.extend(self.extra_topos)
        self.parameters.extend(self.extra_params)

        # XXX - matcher currently not used because I don't trust OPLS AA/M
        # Initialize matcher only now that all topologies etc have been set
        if self.forcefield in []:
            self.matcher = GromacsMatcher(topologies=self.topologies)

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
                                  hmr=self.hmr,
                                  extra_topos=self.extra_topos,
                                  extra_params=self.extra_params,
                                  override_defaults=self.override)
            print("Writing intermediate psf")
            psfgen.write(self.outprefix)
            self._psf_to_gromacs()

        elif "amber" in self.forcefield:
            prmgen = AmberWriter(molid=self.molid,
                                 tmp_dir=self.tmp_dir,
                                 forcefield=self.forcefield,
                                 hmr=self.hmr,
                                 lipid_sel=self.lipid_sel,
                                 extra_topos=self.extra_topos,
                                 extra_params=self.extra_params,
                                 override_defaults=self.override)
            print("Writing intermediate prmtop")
            prmgen.write(self.outprefix)
            self._amber_to_gromacs()

        # Now native GROMACS style for gromos or opls
        else:
            # Currently unsupported
            raise DabbleError("Forcefield '%s' not supported for gromacs"
                              % self.forcefield)

            print("Using the following topology files and/or directories:")
            for top in self.topologies:
                print("  - %s" % os.path.split(top)[1])

            self._set_atom_names()
            self._run_pdb2gmx()

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

    def _set_atom_names(self):
        """
        Sets the atom and residue names for a GROMACS-style force field

        Args:
            molid (int): Molecule ID
        """
        # Rename all residues
        residues = set(atomsel("all", molid=self.molid).residue)
        n_res = len(residues)

        while residues:
            if len(residues) % 500 == 0:
                print("Renaming residues.... %.0f%%  \r"
                      % (100.-100*len(residues)/float(n_res)), flush=True)

            residue = residues.pop()
            sel = atomsel("residue %s" % residue)
            resnames, atomnames = self.matcher.get_names(sel,
                                                        print_warning=False)
            if not resnames:
                rgraph = self.matcher.parse_vmd_graph(sel)[0]
                self.matcher.write_dot(rgraph, "rgraph.dot")
                raise ValueError("Unknown residue %s:%d"
                                 % (sel.resname[0], sel.resid[0]))

            self._apply_naming_dictionary(resnames=resnames,
                                          atomnames=atomnames)
            sel.user = 1.0

        # TODO: handle disulfides, other noncanonical bonds

    #==========================================================================

    def _run_pdb2gmx(self):
        """
        Runs pdb2gmx, creating a topology and structure file

        Args:
            TODO

        Returns:
            (str) Prefix of file written

        Raises:
            DabbleError if pdb2gmx is not present
        """

        # Ensure gmx is actually available
        try:
            subprocess.call("gmx", stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            raise DabbleError("gmx executable not found. Is GROMACS installed "
                              "and accessible in $PATH?")

        # Save the system as a .gro file
        _, grout = tempfile.mkstemp(prefix=self.outprefix,
                                    suffix=".gro", dir=self.tmp_dir)
        os.close(_)
        atomsel("all", self.molid).write("gro", grout)

        # Prepare arguments to pdb2gmx
        gmx_args = [
                    "gmx", "pdb2gmx",
                    "-f", grout,
                    "-o", self.outprefix + ".gro",
                    "-p", self.outprefix + ".top",
                    "-ff", os.path.split(self.ffdir)[1],
                    "-water", "tip3p", # TODO water models
                    "-noignh", # Use hydrogens in coordinate file
                    ]

        if self.hmr:
            gmx_args += ["-heavyh"]

        # Set GMXLIB environment variable to the directory containing our
        # forcefield, or the user provided one, so pdb2gmx will find it
        oldlib = os.environ.get("GMXLIB", "")
        os.environ["GMXLIB"] = os.path.split(self.ffdir)[0]

        # Run pdb2gmx
        out = subprocess.check_output(gmx_args).decode("utf-8")
        try:
            out = "\n==============BEGIN PDB2GMX OUTPUT==============\n" + out \
                + "\n===============END PDB2GMX OUTPUT===============\n"

            if self.debug:
                print(out)

        except:
            print(out)
            raise DabbleError("Call to pdb2gmx failed! See above output for "
                              "more information")

        # Clean up
        os.environ["GMXLIB"] = oldlib

    #==========================================================================
    #                              Static methods
    #==========================================================================

    @classmethod
    def get_topologies(cls, forcefield):
        """
        Gets the path to GROMACS-format topologies for a given force field
        """

        # Amber, Charmm, and OPLS handled by conversion
        if forcefield == "charmm":
            return CharmmWriter.get_topologies(forcefield)

        elif forcefield == "amber":
            return AmberWriter.get_topologies(forcefield)

        elif forcefield == "opls":
            return CharmmWriter.get_topologies(forcefield)

        # No forcefields really ship with gromacs right now because
        # I found an error in the OPLS AA/M gromacs implementation and
        # they won't respond to my emails

        # Use GROMACS forcefield for the remaining ones
        #if forcefield == "opls":
        #    ffdir = "oplsaam.ff"

        #elif forcefield == "gromos":
        #    ffdir = "gromos54a7.ff"

        else:
            raise DabbleError("Unsupported forcefield %s" % forcefield)

        #return [cls.get_forcefield_path(ffdir)]

    #==========================================================================

    @classmethod
    def get_parameters(cls, forcefield):
        """
        Get the path to GROMACS-format parameter files for a given force field
        """

        # Amber and Charmm handled by converstion
        if forcefield == "charmm":
            return CharmmWriter.get_parameters(forcefield)

        elif forcefield == "amber":
            return AmberWriter.get_parameters(forcefield)

        elif forcefield == "opls":
            return CharmmWriter.get_parameters(forcefield)

        # Gromacs topologies and parameters are the same directories
        else:
            return GromacsWriter.get_topologies(forcefield)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
