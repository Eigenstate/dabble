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
import tempfile

from abc import ABC, abstractmethod
from dabble import DabbleError
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
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                                CONSTANTS                                    #
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Will be populated by appropriate writer
    WATER_NAMES = {}
    WATER_O_NAME = None
    WATER_H_NAMES = None

    #==========================================================================

    def __init__(self, molid, **kwargs):
        """
        Creates a molecule writer

        Args:
            molid (int): VMD molecule ID of system to write
            tmp_dir (str): Directory for temporary files. Defaults to "."
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            forcefield (str): Force field to use
            water_model (str): Water model to use
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

        self.forcefield = kwargs.get("forcefield")
        self.water_model = kwargs.get("water_model")

        self.extra_topos = kwargs.get("extra_topos")
        self.extra_params = kwargs.get("extra_params")
        # Handle None from argparse in command line invocation
        if self.extra_topos is None:
            self.extra_topos = []
        if self.extra_params is None:
            self.extra_params = []

        # Handle None from argparse in command line invocation
        if self.extra_topos is None:
            self.extra_topos = []
        if self.extra_params is None:
            self.extra_params = []

    #==========================================================================

    def _rename_by_resname(self, resname, renumber=False):
        """
        Tries to match up all of one resname in one go, by getting an
        isomorphism to one residue in the selection and applying the resulting
        name matching to all other residues. Lots of error checking, will
        refuse to rename residues that don't have the same original atom
        names. Only one isomorphism will be  performed, non-matching residues
        with this resname will be left alone.

        Args:
            resname (str): Residue name to rename
            renumber (bool): Whether to renumber the resids of renamed residues

        Returns:
            (list of int): Residues that were successfully matched
        """
        residues = set(atomsel("resname '%s'" % resname,
                               molid=self.molid).residue)
        matched = []

        # Gather info about first residue. We'll check the others against this.
        firstres = residues.pop()
        sel = atomsel('residue %s' % firstres)
        resid = 1

        # Find a name match for the first residue
        newresname, atomnames = self.matcher.get_names(sel, print_warning=True)
        if not newresname:
            raise DabbleError("No residue definition for residue %s:%s, "
                              "residue %d"
                              % (sel.name[0], sel.resid[0], sel.residue[0]))

        refnames = sel.name # Correct names, before renaming

        # Do the renaming for this residue
        self._apply_naming_dictionary(atomnames=atomnames,
                                      resnames=newresname,
                                      verbose=False)
        newnames = sel.name # New names, after renaming
        newresname = sel.resname[0]
        matched.append(firstres)

        if renumber:
            sel.resid = resid
            resid += 1

        # Now loop through all other residues in the selection
        for res in residues:
            sel = atomsel("residue %s" % res)

            # Check original atom names are the same
            if sel.name != refnames:
                logger.warning("Resid %s doesn't appear to be a standard "
                               "'%s' residue. Will attempt to match another "
                               "topology definition." % (sel.resid[0], resname))
                continue

            # Match up atom names. Since we know atoms are in the same
            # order as the reference residue, we can safely assign refatoms
            # to our names without more atom selections
            sel.name = newnames
            sel.resname = newresname
            matched.append(res)

            if renumber:
                sel.resid = resid
                resid += 1

        return matched

    #==========================================================================

    def _set_water_names(self):
        """
        Sets the names of water residues and atoms according to the given
        water model. We do it this way instead of with the GraphMatcher because
        waters can have a fake bond
        """
        # Sanity check
        if self.water_model not in self.WATER_NAMES:
            raise DabbleError("Unsupported water model '%s' with forcefield "
                              "'%s'" % (self.water_model, self.forcefield))

        watres = self.WATER_NAMES[self.water_model]
        if watres not in self.matcher.known_res:
            raise DabbleError("Water resname '%s' for model '%s' not defined "
                              "in topology files" % (watres, self.water_model))

        # Set consistent residue and atom names, crystal waters
        # can be named HOH, etc
        residues = set(atomsel("water").residue)

        # If no water, nothing to do
        if not residues:
            return

        watsel = "residue %s" % ' '.join(str(_) for _ in residues)

        atomsel(watsel).resname = self.WATER_NAMES[self.water_model]
        atomsel("%s and noh" % watsel).name = self.WATER_O_NAME
        atomsel("%s and not noh" % watsel).name = self.WATER_H_NAMES * len(residues)

    #==========================================================================

    def _write_water_pdbs(self):
        """
        Renames waters and writes them, 10,000 at a time, to PDB files.
        This handles character limit in PDB files, and multiple water
        models

        Returns:
            (list of str): PDB files saved
        """

        allw = atomsel('water and user 1.0')
        pdbs = []
        print("Found %d water residues" % len(set(allw.residue)))

        # Find the problem waters with unordered indices
        problems = []
        for r in set(allw.residue):
            widx = atomsel('residue %s' % r).index
            if max(widx) - min(widx) != 2:
                problems.append(r)
                atomsel('residue %s' % r).user = 0.0 # get it out of allw

        allw.update()
        num_written = int(len(allw)/(9999*3))+1
        print("Going to write %d files for %d water atoms"
              % (num_written, len(allw)))

        # Pull out and write 10k waters at a time if we have normal waters
        if allw:
            for i in range(num_written):
                _, temp = tempfile.mkstemp(suffix='_%d.pdb' % i,
                                           prefix='psf_wat_',
                                           dir=self.tmp_dir)
                os.close(_)
                residues = list(set(allw.residue))[:9999]

                batch = atomsel("residue %s"% ' '.join(str(x) for x in residues))
                try:
                    batch.resid = [k for k in range(1, int(len(batch)/3)+1)
                                   for _ in range(3)]
                except ValueError:
                    raise DabbleError("\nERROR! You have some waters missing "
                                      "hydrogens!\nFound %d water residues, but"
                                      " %d water atoms. Check your crystal "
                                      "waters in the input structure."
                                      % (len(residues), len(batch)))
                batch.user = 0.0
                batch.write('pdb', temp)
                pdbs.append(temp)
                allw.update()

        # Now write the problem waters
        if problems:
            updb = self._write_unorderedindex_waters(problems, self.molid)
            pdbs.append(updb)

        return pdbs

    #==========================================================================

    def _write_unorderedindex_waters(self, residues, molid):
        """
        Renumbers and sorts the specified waters manually. This is much less
        efficient but is necessary in cases where atoms within a water molecule
        are not sequential in index, preventing quick renaming with VMD.
        Identify problem waters, then call this on them. It'll write its own
        psf_wat_* file with just those waters, minimizing inefficiency.

        Args:
            residues (list of int): Problem water molecules
            molid (int): VMD molecule ID to write
        Returns:
            (str): Filename where waters are written
        """
        f, temp = tempfile.mkstemp(suffix='_indexed.pdb', prefix='psf_wat_',
                                   dir=self.tmp_dir)
        idx = 1
        with os.fdopen(f, 'w') as fileh:
            for ridx, residue in enumerate(residues):
                res = atomsel('residue %d' % residue, molid=molid)

                for i in res.index:
                    a = atomsel('index %d' % i, molid)
                    fileh.write(self.get_pdb_line(a, idx, ridx+1))
                    idx += 1

            fileh.write('END\n')
        return temp

    #==========================================================================

    @abstractmethod
    def write(self, filename):
        pass

    @abstractmethod
    def get_topologies(self, forcefield):
        pass

    @abstractmethod
    def get_parameters(self, forcefield):
        pass

    #==========================================================================
    #                              Static methods
    #==========================================================================

    @staticmethod
    def _apply_naming_dictionary(resnames, atomnames, verbose=False):
        """
        Applies the atom names from a matcher.

        Args:
            resnames (dict int->str or str): Atom index to residue name,
                or just residue name for all atoms
            atomnames (dict int->str): Atom index to atom name
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
                                     ins if ins else " ",
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
