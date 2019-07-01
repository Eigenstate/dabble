"""
This module contains the CharmmWriter class and associated methods,
which outputs a psf/pdb file with CHARMM names and parameters.
It does this by converting atom names to CHARMM names, writing
intermediate files as necessary to invoke the vmd psfgen plugin.

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
import tempfile
from parmed.formats.registry import load_file
from psfgen import PsfGen
from vmd import atomsel, molecule

from dabble import DabbleError
from dabble.param import CharmmMatcher, MoleculeWriter, Patch

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                CONSTANTS                                    #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

patchable_acids = ('ACE ALA ARG ASN ASP CYS CYX GLN GLU GLY HIE HIS HSP HSE '
                   'HSD ILE LEU LYS MET NMA PHE PRO SER THR TRP TYR VAL')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmWriter(MoleculeWriter):
    """
    An object that handles all the conversions to a psf file
    by interfacing with psfgen.

    Writes a pdb/psf file pair from the current molecule using the
    CHARMM36 topology and atom names/types. Interfaces with psfgen by
    dynamically generating the .tcl file that psfgen takes as input.
    Prompts the user for additional topology files and helps with
    matching atom names that cannot be automatically translated to the
    charmm naming conventions.
    """

    #==========================================================================

    def __init__(self, molid, **kwargs):
        """
        Creates a CHARMM writer

        Args:
            molid (int): VMD molecule ID of system to write
            tmp_dir (str): Directory for temporary files. Defaults to "."
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            hmr (bool): If hydrogen masses should be repartitioned. Defaults
                to False
            forcefield (str): Forcefield to use, either "charmm" or "amber"
            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default forcefield parameters.
            debug_verbose (bool): Prints additional output, like from psfgen.
        """

        # Initialize default options
        super(CharmmWriter, self).__init__(molid, **kwargs)

        # Create a psf generator object
        self.psfgen = PsfGen()

        # Set forcefield default topologies and parameters
        self.forcefield = kwargs.get("forcefield", "charmm")
        self.topologies = self.get_topologies(self.forcefield)
        self.parameters = self.get_parameters(self.forcefield)

        if "charmm" in self.forcefield:
            if self.hmr:
                raise DabbleError("HMR not supported with CHARMM ff yet")

        # Handle override and extra topologies
        if self.override:
            self.topologies = []
            self.parameters = []

        # Now extra topologies (put in self by super __init__)
        self.topologies.extend(self.extra_topos)
        self.parameters.extend(self.extra_params)

        # Once all topologies defined, initialize matcher only if
        # using CHARMM topologies (not if we're doing a conversion)
        if "charmm" in self.forcefield or "opls" in self.forcefield:
            self.matcher = CharmmMatcher(self.topologies)

    #=========================================================================

    def write(self, filename):
        """
        Writes the parameter and topology files

        Args:
            filename (str): File name to write. File type suffix will be added.
        """
        self.outprefix = filename

        # Put our molecule on top
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        # Amber forcefield done with AmberWriter then conversion
        if "amber" in self.forcefield:
            # Avoid circular import by doing it here
            from dabble.param import AmberWriter
            prmtopgen = AmberWriter(molid=self.molid,
                                    tmp_dir=self.tmp_dir,
                                    forcefield=self.forcefield,
                                    hmr=self.hmr,
                                    lipid_sel=self.lipid_sel,
                                    extra_topos=self.extra_topos,
                                    extra_params=self.extra_params,
                                    override_defaults=self.override,
                                    debug_verbose=self.debug)
            prmtopgen.write(self.outprefix)
            self._prmtop_to_charmm()

        # Charmm forcefield
        elif "charmm" in self.forcefield:
            self._run_psfgen()

        # OPLS forcefield. Same as charmm but list separately for readability
        elif "opls" in self.forcefield:
            self._run_psfgen()

        else:
            raise DabbleError("Unsupported forcefield '%s' for CharmmWriter"
                              % self.forcefield)

        # Check output and finish up
        self._check_psf_output()

        # Reset top molecule
        molecule.set_top(old_top)


    #=========================================================================
    #                           Static methods                               #
    #=========================================================================

    @classmethod
    def get_topologies(cls, forcefield):

        if forcefield == "charmm":
            topos = [
                "top_all36_caps.rtf",
                "top_all36_cgenff.rtf",
                "top_all36_prot.rtf",
                "top_all36_lipid.rtf",
                "top_all36_carb.rtf",
                "top_all36_na.rtf",
                "toppar_water_ions.str",
                "toppar_all36_prot_na_combined.str",
                "toppar_all36_prot_fluoro_alkanes.str"
            ]

        elif forcefield == "opls":
            topos = [
                "opls_aam.rtf",
                "opls_aam_caps.rtf"
            ]

        elif forcefield == "amber":
            from dabble.param import AmberWriter # avoid circular dependency
            return AmberWriter.get_topologies(forcefield)

        else:
            raise ValueError("Invalid forcefield: '%s'" % forcefield)

        return [cls._get_forcefield_path(top) for top in topos]

    #=========================================================================

    @classmethod
    def get_parameters(cls, forcefield):

        if forcefield == "charmm":
            prms = [
                "toppar_water_ions.str",
                "par_all36m_prot.prm",
                "par_all36_cgenff.prm",
                "par_all36_lipid.prm",
                "par_all36_carb.prm",
                "par_all36_na.prm",
                "toppar_all36_prot_na_combined.str"
            ]

        elif forcefield == "amber":
            from dabble.param import AmberWriter # avoid circular dependency
            return AmberWriter.get_parameters(forcefield)

        elif forcefield == "opls":
            prms = [
                "opls_aam.prm"
            ]

        else:
            raise ValueError("Invalid forcefield: '%s'" % forcefield)

        return [cls._get_forcefield_path(par) for par in prms]

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _write_water_blocks(self):
        """
        Writes a lot of temporary files with 10000 waters each, to bypass
        psfgen being stupid with files containing more than 10000 of a residue.
        """
        # Put current molecule on top to simplify atom selection language
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        # Set consistent residue and atom names, crystal waters
        # can be named HOH, etc
        atomsel('water').resname = "TIP3"
        atomsel('resname TIP3').chain = "W"
        atomsel('resname TIP3 and element O').name = "OH2"

        # Dowser can name water hydrogens strangely
        atomsel('resname TIP3 and name HW1').name = "H1"
        atomsel('resname TIP3 and name HW2').name = "H2"

        # Select all the waters. We'll use the user field to track which
        # ones have been written
        allw = atomsel('water and user 1.0')
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
                temp = tempfile.mkstemp(suffix='_%d.pdb' % i, prefix='psf_wat_',
                                        dir=self.tmp_dir)[1]
                residues = list(set(allw.residue))[:9999]

                batch = atomsel('residue %s' % ' '.join([str(x) for x in residues]))
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
                allw.update()

                self.psfgen.add_segment(segid="W%d" % i, pdbfile=temp)
                self.psfgen.read_coords(segid="W%d" % i, filename=temp)

        # Now write the problem waters
        updb = self._write_unorderedindex_waters(problems, self.molid)
        self.psfgen.add_segment(segid="W%d" % (num_written+1), pdbfile=updb)
        self.psfgen.read_coords(segid="W%d" % (num_written+1), filename=updb)

        molecule.set_top(old_top)

        return num_written

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
        temp = tempfile.mkstemp(suffix='_indexed.pdb', prefix='psf_wat_',
                                dir=self.tmp_dir)[1]
        fileh = open(temp, 'w')

        idx = 1
        for ridx, residue in enumerate(residues):
            res = atomsel('residue %d' % residue, molid=molid)

            for i in res.index:
                a = atomsel('index %d' % i, molid) # pylint: disable=invalid-name
                fileh.write(self.get_pdb_line(a, idx, ridx+1))
                idx += 1

        fileh.write('END\n')
        fileh.close()
        return temp

    #==========================================================================

    def _write_lipid_blocks(self):
        """
        Writes a temporary PDB file containing the lipids for later use by
        psfgen. Renumbers the lipid residues because some can have **** instead
        of an integer for resid in large systems, which will crash psfgen. Also
        sets atom names for some common lipids (currently POPC)

        Raises:
            NotImplementedError if more than 10,000 lipids are present since it
              doesn't support feeding multiple lipid blocks to psfgen currently
            NotImplementedError if lipid other than POPC,POPE,POPG is found
        """
        # Put current molecule on top to simplify atom selection
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        # Collect lipid residues up
        alll = atomsel('(%s) and user 1.0' % self.lipid_sel)
        residues = list(set(alll.residue))
        residues.sort()

        # Lipids not compatible with AMBER parameters, CHARMM format
        if len(residues) and \
                ("amber" in self.forcefield or "opls" in self.forcefield):
            raise ValueError("AMBER or OPLS parameters not supported for lipids"
                             " in CHARMM output format")

        # Sanity check for < 10k lipids
        if len(residues) >= 10000:
            raise NotImplementedError("More than 10k lipids found")

        # Loop through all residues and renumber and correctly name them
        counter = 1
        for res in residues:
            sel = atomsel('residue %s' % res)

            # Find a name match for this residue
            (newname, atomnames) = self.matcher.get_names(sel, print_warning=True)
            if not newname:
                raise DabbleError("No residue definition for lipid %s:%s"
                                  % (sel.name[0], sel.resid[0]))

            # Do the renaming
            self._apply_naming_dictionary(atomnames, resnames=newname,
                                          verbose=True)
            # Renumber residue
            sel.resid = counter
            counter = counter + 1


        # Write temporary lipid pdb
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_lipid_',
                                dir=self.tmp_dir)[1]
        alll.user = 0.0
        alll.write('pdb', temp)

        # Generate lipid segment
        self.psfgen.add_segment(segid="L", pdbfile=temp)
        self.psfgen.read_coords(segid="L", filename=temp)

        # Put old top back
        molecule.set_top(old_top)

    #==========================================================================

    def _write_ion_blocks(self):
        """
        Writes a PDB file containing correctly named ions for use by
        psfgen, and instructs psfgen to use it in TCL code.
        """

        # Put our molecule on top to simplify atom selection language
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        # Get ion resids that aren't associated w other molecules
        # because some ligands have Na, Cl, K
        total = atomsel('element Na Cl K')
        if len(total):
            not_ions = atomsel("(same fragment as element Na Cl K)  and (not index %s)"
                               % " ".join([str(s) for s in set(total.index)]))
            ions = set(total.residue) - set(not_ions.residue)
        else:
            ions = set(total.residue)

        if not len(ions):
            return
        ionstr = "residue " + " ".join([str(s) for s in ions])

        # Fix the names
        atomsel('%s and name NA' % ionstr).name = "SOD"
        atomsel('%s and name CL' % ionstr).name = "CLA"
        atomsel('%s and name K' % ionstr).name = "POT"
        atomsel('%s and name NA' % ionstr).resname = "SOD"
        atomsel('%s and name CL' % ionstr).resname = "CLA"
        atomsel('%s and name K' % ionstr).resname = "POT"

        # Renumber the residues since some may be above 10k
        residues = atomsel('name SOD CLA POT').residue
        batch = atomsel('residue %s' % ' '.join([str(s) for s in set(residues)]))
        batch.resid = [k for k in range(1, len(batch)+1)]

        # Save the temporary ions file
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_ions_',
                                dir=self.tmp_dir)[1]
        atomsel('name SOD CLA POT').user = 0.0
        atomsel('name SOD CLA POT').write('pdb', temp)

        self.psfgen.add_segment(segid="I", pdbfile=temp)
        self.psfgen.read_coords(segid="I", filename=temp)

        molecule.set_top(old_top)

    #==========================================================================

    def _find_single_residue_names(self, resname, molid):
        """
        Uses graph matcher and available topologies to match up
        ligand names automatically. Tries to use graphs, and if there's an
        uneven number of atoms tries to match manually to suggest which atoms
        are most likely missing.

        Args:
          resname (str): Residue name of the ligand that will be written.
            All ligands will be checked separately against the graphs.
          molid (int): VMD molecule ID to consider

        Returns:
          (list of ints): Residue numbers (not resid) of all input ligands
            that were successfully matched. Need to do it this way since
            residue names can be changed in here to different things.

        Raises:
          ValueError if number of resids does not match number of residues as
            interpreted by VMD
          NotImplementedError if a residue could not be matched to a graph.
        """
        # Put our molecule on top
        old_top = molecule.get_top()
        molecule.set_top(molid)

        # Sanity check that there is no discrepancy between defined resids and
        # residues as interpreted by VMD.
        for chain in set(atomsel("user 1.0 and resname '%s'" % resname).chain):
            residues = list(set(atomsel("user 1.0 and resname '%s' and chain %s"
                                        % (resname, chain)).residue))
            resids = list(set(atomsel("user 1.0 and resname '%s' and chain %s"
                                      % (resname, chain)).resid))
            if len(residues) != len(resids):
                raise DabbleError("VMD found %d residues for resname '%s', "
                                  "but there are %d resids! Check input."
                                  % (len(residues), resname,
                                     len(resids)))

        for residue in residues:
            sel = atomsel("residue %s and resname '%s' and user 1.0" % (residue, resname))

            (newname, atomnames) = self.matcher.get_names(sel, print_warning=True)
            if not newname:
                (resname, patch, atomnames) = self.matcher.get_patches(sel)

                if not newname:
                    print("ERROR: Could not find a residue definition for %s:%s"
                          % (resname, residue))
                    raise NotImplementedError("No residue definition for %s:%s"
                                              % (resname, residue))
                print("\tApplying patch %s to ligand %s" % (patch, newname))

            # Do the renaming
            self._apply_naming_dictionary(atomnames, resnames=newname,
                                          verbose=True)

        #logger.info("Renamed %d atoms for all resname %s->%s" % (num_renamed, resname, name))
        molecule.set_top(old_top)

        return residues

    #==========================================================================

    def _write_generic_block(self, residues):
        """
        Matches ligands to available topology file, renames atoms, and then
        writes temporary files for the ligands

        Args:
          residues (list of int): Residue numbers to be written. Will all
            be written to one segment.

        Returns:
          True if successful
        """
        # Put our molecule on top to simplify atom selection language
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        alig = atomsel('user 1.0 and residue %s' % " ".join([str(x) for x in residues]))

        # Write temporary file containg the residues and update tcl commands
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_block_',
                                dir=self.tmp_dir)[1]
        alig.write('pdb', temp)
        alig.user = 0.0

        self.psfgen.add_segment(segid="B%s" % residues[0], pdbfile=temp)
        self.psfgen.read_coords(segid="B%s" % residues[0],
                                filename=temp)

        if old_top != -1:
            molecule.set_top(old_top)
        return True

    #==========================================================================

    def _write_protein_blocks(self, molid, frag):
        """
        Writes a protein fragment to a pdb file for input to psfgen
        Automatically assigns amino acid names

        Args:
            molid (int): VMD molecule ID of renumbered protein
            frag (str): Fragment to write

        Returns:
            (list of Patches): Patches to add to psfgen input files
       """

        print("Setting protein atom names")

        # Put our molecule on top to simplify atom selection language
        old_top = molecule.get_top()
        molecule.set_top(molid)
        patches = set()
        extpatches = set()
        seg = "P%s" % frag

        residues = list(set(atomsel("fragment '%s'" % frag).residue))
        for residue in residues:
            sel = atomsel('residue %s' % residue)
            resid = sel.resid[0]

            # Only try to match single amino acid if there are 1 or 2 bonds
            if len(self.matcher.get_extraresidue_atoms(sel)) < 3:
                (newname, atomnames) = self.matcher.get_names(sel, False)

            # See if it's a disulfide bond participant
            else:
                (newname, patch, atomnames) = \
                        self.matcher.get_disulfide("residue %d" % residue,
                                                   frag, molid)
                if newname:
                    extpatches.add(patch)

            # Couldn't find a match. See if it's a patched residue
            if not newname:
                (newname, patchname, atomnames) = self.matcher.get_patches(sel)
                if newname:
                    # This returns patch name only, not a Patch object
                    patches.add(Patch(name=patchname, segids=[seg],
                                      resids=[resid]))

            # Fall through to error condition
            if not newname:
                raise DabbleError("Couldn't find a patch for %s:%s"
                                  % (sel.resname[0], resid))

            # Do the renaming
            self._apply_naming_dictionary(atomnames, resnames=newname)

        # Save protein chain in the correct order
        filename = self.tmp_dir + '/psf_protein_%s.pdb' % seg
        _write_ordered_pdb(filename, "fragment '%s'" % frag, molid)
        print("\tWrote %d atoms to the protein segment %s"
              % (len(atomsel("fragment %s" % frag)), seg))

        # Now invoke psfgen for the protein segments
        self.psfgen.add_segment(segid=seg, pdbfile=filename)

        print("Applying the following single-residue patches to P%s:\n" % frag)
        print("\t%s" % "\t".join(str(_) for _ in patches))
        for p in patches:
            self.psfgen.patch(patchname=p.name, targets=p.targets())

        self.psfgen.read_coords(segid=seg, filename=filename)

        if old_top != -1:
            molecule.set_top(old_top)

        return extpatches

    #==========================================================================

    def _check_psf_output(self):
        """
        Scans the output psf from psfgen for atoms where the coordinate
        could not be set, indicating an unmatched atom. This check is necessary
        because sometimes psfgen will run with no errors or warnings but will
        have unmatched atoms that are all at (0,0,0).
        """

        # Check file was written at all
        if not os.path.isfile('%s.pdb'% self.outprefix):
            raise DabbleError("\nERROR: psf file failed to write.\n"
                              "       Please see log above.\n")

        # Open the pdb file in VMD and check for atoms with no occupancy
        fileh = molecule.load('pdb', '%s.pdb' % self.outprefix)
        errors = atomsel("occupancy=-1", molid=fileh)

        # Print out error messages
        if errors:
            errstr = "\nERROR: Couldn't find the following atoms."
            for i in range(len(errors)):
                errstr += "  %s%s:%s" % (errors.resname[i], errors.resid[i],
                                         errors.name[i])

            errstr += "Check if they are present in the original structure.\n"
            raise DabbleError(errstr)
        else:
            print("\nChecked output pdb/psf has all atoms present "
                  "and correct.\n")


    #==========================================================================

    def _find_residue_in_rtf(self, resname, molid):
        """
        Scans the input topology files to find a name match for the given
        residue name, then pulls out the atoms involved and checks that they
        are all present in the input coordinates, prompting the user to correct
        the names of atoms that could not be matched.

        Residue ID is used because there can be multiple copies of a residue
        with the same name, but only one has missing or extra atoms.

        Args:
          resname (str): Residue name to check
          molid (int): VMD molecule ID

        Returns:
          True if all matching was successful
          False if the residue name cannot be found
        """

        print("Finding residue name '%s'" % resname)
        for top in self.topologies:
            topfile = open(top, 'r')
            topo_atoms = _get_atoms_from_rtf(text=topfile.readlines(),
                                             resname=resname)
            # Use first definition found of this residue
            if len(topo_atoms):
                break
            topfile.close()
        if not len(topo_atoms):
            return False
        print("Successfully found residue %s in input topologies" % resname)

        # Match up atoms with python sets
        pdb_atoms = set(atomsel("resname '%s' and user 1.0" % resname,
                                molid=molid).name)
        pdb_only = pdb_atoms - topo_atoms
        topo_only = topo_atoms - pdb_atoms

        # If uneven number of atoms, there are missing or additional atoms
        if len(pdb_atoms) > len(topo_atoms):
            raise DabbleError("\nERROR: Cannot process modified residue %s.\n"
                              "There are %d extra atoms in the input structure "
                              "that are undefined in the topology file. The "
                              "following atoms could not be matched and may "
                              "either be misnamed, or additional atoms:\n"
                              "[ %s ]\n"
                              % (resname, len(pdb_atoms)-len(topo_atoms),
                                 " ".join(pdb_only)))

        if len(topo_atoms) > len(pdb_atoms):
            raise DabbleError("\nERROR: Cannot process modified residue %s.\n"
                              "There are %d missing atoms in the input structure "
                              "that are defined in the topology file. The "
                              "following atoms could not be matched and may "
                              "either be misnamed or deleted atoms:\n"
                              "[ %s ]\n"
                              % (resname, len(topo_atoms)-len(pdb_atoms),
                                 " ".join(topo_only)))

        # Offer to rename atoms that couldn't be matched to the topology
        if len(pdb_only):
            print("\nWARNING: Having some trouble with modified residue %s.\n"
                  "         The following atom names cannot be matched up "
                  " to the input topologies. They are probably "
                  " misnamed.\n" % resname)
            print("         To help you, here are the atom names that "
                  " should be present according to the topology "
                  " but were not found:\n")
            print("         [ %s ]\n" % ' '.join([str(t) for t in topo_only]))
            print(" Please enter a valid name for each atom as "
                  "it appears or CTRL+D to quit..\n")
            for unmatched in pdb_only:
                print("Unmatched topology names: [ %s ]"
                      % ' '.join(topo_only))

                try: input = raw_input
                except NameError: pass
                newname = input("  %s  -> " % unmatched)
                while newname not in topo_only:
                    print("'%s' is not an available name in the topology."
                          "Please try again.\n" % newname)
                    newname = input("  %s  -> " % unmatched)
                atomsel("resname '%s' and user 1.0 and name '%s'"
                        % (resname, unmatched)).name = newname
                pdb_atoms = set(atomsel("resname '%s' and user 1.0"
                                        % resname).name)
                topo_only = topo_atoms-pdb_atoms
                resname = newname

            # Recurse to check that everything is assigned correctly
            self._find_residue_in_rtf(resname, molid)
        print("Matched up all atom names for resname '%s'\n" % resname)
        return True

    #==========================================================================

    def _get_patch(self, seg, resid):
        """
        Prompts the user for a patch to apply for the given residue.
        Gathers available patches from topology files

        Args:
          seg (str): Segment to apply the patch to
          resid (int): Residue ID to apply the patch to

        Returns:
          (str) patch line to put in the psfgen input file
        """
        avail_patches = self._get_avail_patches()
        print("What is the patch name I should apply?")
        print("Type NONE for no patch, if your residue is completely "
              "defined in a str file")
        print("Or type HELP for a list of all patches I know about")
        try: input = raw_input
        except NameError: pass
        patchname = input("> ")
        if patchname == "HELP":
            print("   PATCH     COMMENT")
            print("   -----     -------")
            for patch in avail_patches:
                print("%7s %s" % (patch, avail_patches[patch]))
            patchname = input("> ")
        while (patchname not in avail_patches) and (patchname != "NONE"):
            print("I don't know about patch %s" % patchname)
            patchname = input("Try again > ")
        if patchname == "NONE":
            return ""

        return "patch %s %s:%d\n" % (patchname, seg, resid)

    #==========================================================================

    def _get_avail_patches(self):
        """
        Gathers the patches defined in all topology files.

        Returns:
          (dict str -> str): Patch names as keys, comment as value
        """
        avail_patches = {}
        for top in self.topologies:
            topfile = open(top, 'r')
            for line in topfile:
                tokens = line.split()
                if not len(tokens):
                    continue
                if tokens[0] == "PRES":
                    comment = ' '.join(tokens[tokens.index("!")+1:])
                    avail_patches[tokens[1]] = comment
        return avail_patches

    #==========================================================================

    def _run_psfgen(self):

        # Read topology files in to psfgen
        print("Using the following topologies:")
        for top in self.topologies:
            print("  - %s" % os.path.split(top)[1])
            self.psfgen.read_topology(top)


        # Mark all atoms as unsaved with the user field
        atomsel('all', molid=self.molid).user = 1.0
        check_atom_names(molid=self.molid)

        # Now ions if present, changing the atom names
        if len(atomsel('ions', molid=self.molid)) > 0:
            self._write_ion_blocks()

        # Save water 10k molecules at a time
        if len(atomsel('water', molid=self.molid)):
            self._write_water_blocks()

        # Now lipid
        # TODO: lipid names aren't matched with topologies?
        if len(atomsel(self.lipid_sel)):
            self._write_lipid_blocks()

        # Now handle the protein
        # Save and reload the protein so residue looping is correct
        if len(atomsel("resname %s" % patchable_acids, molid=self.molid)):
            extpatches = set()
            for frag in sorted(set(atomsel("resname %s" % patchable_acids,
                                    molid=self.molid).fragment)):
                extpatches.update(self._write_protein_blocks(self.molid, frag))
            atomsel("same fragment as resname %s" % patchable_acids,
                    molid=self.molid).user = 0.0

            # List all patches applied to the protein
            print("Applying the following patches:\n")
            print("\t%s" % "\n\t".join(str(_) for _ in extpatches))

            # Apply all multi segment patches to the protein
            for p in extpatches:
                self.psfgen.patch(p.name, p.targets())
        else:
            print("\tDidn't find any protein. Continuing...\n")

        # Regenerate angles and dihedrals after applying patches
        # Angles must be regenerated FIRST!
        # See http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2009-2010/4137.html
        self.psfgen.regenerate_angles()
        self.psfgen.regenerate_dihedrals()

        # Check if there is anything else and let the user know about it
        leftovers = atomsel('user 1.0', molid=self.molid)
        for lig in set(leftovers.resname):
            residues = self._find_single_residue_names(resname=lig,
                                                       molid=self.molid)
            self._write_generic_block(residues)

        # Write the output files and run
        self.psfgen.write_psf(filename="%s.psf" % self.outprefix, type="x-plor")
        self.psfgen.write_pdb(filename="%s.pdb" % self.outprefix)

    #==========================================================================

    def _prmtop_to_charmm(self):
        """
        Converts an AMBER prmtop with AMBER parameters to a psf file,
        using ParmEd.
        """
        # Save PSF topology and parameter file
        parmstruct = load_file(self.outprefix + ".prmtop",
                               xyz=self.outprefix + ".inpcrd",
                               structure=True)
        parmstruct.save(self.outprefix + ".psf", format="psf")

        # Save PDB file with coordinates
        m = molecule.load("parm7", self.outprefix + ".prmtop",
                          "rst7", self.outprefix + ".inpcrd")
        atomsel("all", m).write("pdb", self.outprefix + ".pdb")
        molecule.delete(m)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                 FUNCTIONS                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _write_ordered_pdb(filename, sel, molid):
    """
    Writes a pdb file in order of residues, renumbering the atoms
    accordingly, since psfgen wants each residue sequentially while
    VMD will write them in the same order as input, which from Maestro
    created files has some guessed atoms at the end.

    Args:
      filename (str): Name of the pdb file to write
      sel (str): VMD atomsel string for atoms that will be written
      molid (int): VMD molecule ID to write from
    """
    old_top = molecule.get_top()
    molecule.set_top(molid)

    fileh = open(filename, 'w')
    # Use resids since order can be wrong when sorting by residue
    # Then, use residue to pull out each one since it is much much
    # faster then trying to pull out residues
    resids = set(atomsel(sel).resid)

    # Add additional residue constraint to selection since pulling out
    # by resid can match something in a different chain
    resstr = ' '.join([str(x) for x in set(atomsel(sel).residue)])

    idx = 1
    # For renumbering capping groups
    for resid in sorted(resids):
        # Check for alternate locations
        residues = sorted(set(atomsel("resid '%s' and residue %s"
                                      % (resid, resstr)).residue))
        for rid in residues:
            for i in atomsel('residue %d' % rid).index:
                a = atomsel('index %d' % i) # pylint: disable=invalid-name
                fileh.write(MoleculeWriter.get_pdb_line(a, idx, a.resid[0]))
                idx += 1
    fileh.write('END\n')
    atomsel(sel).user = 0.0
    fileh.close()
    molecule.set_top(old_top)

#==========================================================================

def _get_atoms_from_rtf(text, resname):
    """
    Scans the input text for the residue with a given name. Once found,
    pulls out all the atom names that comprise that residue.

    Args:
      text (str): Contents of an rtf file to scan
      resname (str): Residue to look for

    Returns:
      atoms (set of str): Atom names in this residue, or the empyty set if
          the residue was not found.
    """
    atoms = []
    found = False
    for line in text:
        words = line.split()
        if not len(words):
            continue
        if not found and words[0] == 'RESI' \
           and words[1] == resname:
            found = True
        elif found and words[0] == 'ATOM':
            atoms.append(words[1])
        elif found and words[0] == 'RESI':
            break
    return set(atoms)

#==========================================================================

def get_bonded_atoms(molid, index):
    """
    Returns the element of all atoms bonded to the current atom.

    Args:
       molid (int): VMD molecule ID to consider
       index (int): Atom index to look at bonded atoms

    Returns:
      (list of str) elements of atoms bound to the current atom
    """

    asel = atomsel('index %d' % index, molid=molid)
    bound = []
    for atom in asel.bonds[0]:
        bound.append(atomsel('index %d' % atom).element[0])
    return bound

#==========================================================================

def check_atom_names(molid):
    """
    Checks that there are no spaces in atom names. If spaces are
    found, they are removed and a warning is printed
    """

    names = set(atomsel(molid=molid).name)
    for name in names:
        if ' ' in name:
            print("\nWARNING: Found space character in name '%s'\n"
                  "         Incompatible with charmm formats, removing it"
                  % name)
            atomsel("name '%s'", molid=molid).name = name.replace(' ', '')

#==========================================================================
