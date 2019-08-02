"""
 This module contains the AmberWriter class and associated methods,
 which outputs a prmtop/inpcrd file with CHARMM names and parameters
 for use with simulation with the AMBER molecular dynamics package.
 It does this by using the chamber functionality of ParmEd API.

 Author: Robin Betz

 Copyright (C) 2019 Robin Betz
"""
#
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
import warnings
from subprocess import check_output
from vmd import molecule, atomsel

from parmed.tools import chamber, parmout, HMassRepartition, checkValidity
from parmed.amber import AmberParm
from parmed.exceptions import ParameterWarning
from dabble import DabbleError
from dabble.param import MoleculeWriter, AmberMatcher

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberWriter(MoleculeWriter):
    """
    Creates AMBER-format input files, using either CHARMM or AMBER parameters.

    When using the CHARMM parameters, creates a psf/pdb file and interfaces
    with ParmEd chamber command to create AMBER format files.
    When using the AMBER parameters, creates a pdb file and runs the amber
    lipid conversion script to create leap input files.
    """

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #                               CONSTANTS                                  #
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Leaprc handles water name translation for us
    WATER_NAMES = {
                   "tip3"  : "TP3",
                   "tip4e" : "T4E",
                   "spce"  : "SPC",
                  }
    WATER_O_NAME  = "O"
    WATER_H_NAMES = ["H1", "H2"]

    #==========================================================================

    def __init__(self, molid, **kwargs):
        """
        Creates an AMBER Writer.

        Args:
            molid (int): VMD molecule ID of system to write
            tmp_dir (str): Directory for temporary files. Defaults to "."
            forcefield (str): charmm, or amber
            water_model (str): Water model to use
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            hmr (bool): If hydrogen masses should be repartitioned. Defaults
                to False.
            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default amber ff14 parameters.
            debug_verbose (bool): Prints additional output, like from tleap.
        """
        # Initialize standard MoleculeWriter items
        super(AmberWriter, self).__init__(molid, **kwargs)
        self.prompt_params = False

        # Set forcefield default topologies and parameters
        self.forcefield = kwargs.get("forcefield", "amber")
        self.water_model = kwargs.get("water_model", "TIP3")

        if self.forcefield == "charmm":
            # Import CharmmWriter as needed to avoid circular imports
            from dabble.param import CharmmWriter

            # Topologies used will be found and returned by CharmmWriter
            self.topologies = CharmmWriter.get_topologies(self.forcefield,
                                                          self.water_model)
            self.parameters = CharmmWriter.get_parameters(self.forcefield,
                                                          self.water_model)

        # Duplicate code here for clarity
        elif self.forcefield == "opls":
            from dabble.param import CharmmWriter
            self.topologies = CharmmWriter.get_topologies(self.forcefield,
                                                          self.water_model)
            self.parameters = CharmmWriter.get_parameters(self.forcefield,
                                                          self.water_model)

        elif self.forcefield == "amber":
            # This will raise DabbleError if AMBERHOME is unset
            self.topologies = self.get_topologies(self.forcefield,
                                                  self.water_model)
            self.parameters = self.get_parameters(self.forcefield,
                                                  self.water_model)
            self.parameters = []

        else:
            raise DabbleError("Unsupported forcefield: %s" % self.forcefield)

        # Handle override
        if self.override:
            self.topologies = []
            self.parameters = []

        # Now extra topologies and parameters
        self.topologies.extend(self.extra_topos)
        self.parameters.extend(self.extra_params)

        # Initialize matcher only now that all topologies have been set
        if self.forcefield == "amber":
            self.matcher = AmberMatcher(topologies=self.topologies)

    #==========================================================================

    def write(self, filename):
        """
        Writes the parameter and topology files

        Args:
            filename (str): File name to write. File type suffix will be added.
        """
        self.outprefix = filename

        # Charmm forcefield
        if "charmm" in self.forcefield or "opls" in self.forcefield:
            from dabble.param import CharmmWriter
            charmmw = CharmmWriter(molid=self.molid,
                                   tmp_dir=self.tmp_dir,
                                   forcefield=self.forcefield,
                                   water_model=self.water_model,
                                   lipid_sel=self.lipid_sel,
                                   extra_topos=self.extra_topos,
                                   extra_params=self.extra_params,
                                   override_defaults=self.override,
                                   debug_verbose=self.debug)
            charmmw.write(self.outprefix)
            self._psf_to_prmtop()

        # Amber forcefield
        elif "amber" in self.forcefield:

            # Print info about topologies
            print("Using the following topologies:")
            for top in self.topologies:
                print("  - %s" % os.path.split(top)[1])

            # Assign atom types
            print("Assigning AMBER atom types...")
            conect = self._rename_atoms_amber()

            # Create temporary pdb files that will be leap inputs
            pdbs = []
            pdbs.append(self._write_lipids())
            prot_pdbseqs = self._write_protein()
            pdbs.extend(self._write_solvent())
            ligfiles = self._write_ligands()

            # Now invoke leap to create the prmtop and inpcrd
            outfile = self._run_leap(ligfiles, prot_pdbseqs, pdbs,
                                     conect)

            # Repartion hydrogen masses if requested
            if self.hmr:
                print("\nRepartitioning hydrogen masses...")
                parm = AmberParm(prm_name=outfile+".prmtop",
                                 xyz=outfile+".inpcrd")
                action = HMassRepartition(parm, "dowater")
                action.execute()
                write = parmout(action.parm, "%s.prmtop" % self.outprefix)
                write.execute()
                parm = write.parm

            # Check validity of output prmtop using parmed
            parm = AmberParm(prm_name=self.outprefix + ".prmtop",
                             xyz=self.outprefix + ".inpcrd")
            print("\nChecking for problems with the prmtop...")
            print("        Verify all warnings!")
            action = checkValidity(parm)
            action.execute()

        else:
            raise DabbleError("Unhandled forcefield: %s" % self.forcefield)

    #========================================================================#
    #                           Private methods                              #
    #========================================================================#

    def _rename_atoms_amber(self):
        """
        Matches up atom names with those in the provided topologies and
        sets the atom and residue names correctly in the built molecule.
        Handles all non-lipid atoms. Sets the user field of all atoms to 1.0
        to track which things have been written.

        Returns:
            (set of tuples (int,int)): Residue #s of disulfide or otherwise
                noncanonically linked residues

        Raises:
            ValueError if a residue definition could not be found
        """

        self._set_water_names()
        nonlips = set(atomsel("not (water or %s)" % self.lipid_sel,
                              molid=self.molid).residue)
        n_res = len(nonlips)
        conect = set() # Atom indices bound to noncanonical residues
        while nonlips:
            if len(nonlips) % 500 == 0:
                print("Renaming residues.... %.0f%%  \r"
                      % (100.-100*len(nonlips)/float(n_res)), flush=True)

            residue = nonlips.pop()
            sel = atomsel("residue %s" % residue)
            resnames, atomnames = self.matcher.get_names(sel,
                                                         print_warning=False)

            # Check if it's a linkage to another amino acid
            if not resnames:
                resnames, atomnames, other = self.matcher.get_linkage(sel,
                                                                      self.molid)
                if not resnames:
                    rgraph = self.matcher.parse_vmd_graph(sel)[0]
                    self.matcher.write_dot(rgraph, "rgraph.dot")
                    raise DabbleError("ERROR: Could not find a residue definition "
                                      "for %s:%s" % (sel.resname[0],
                                                     sel.resid[0]))

                print("\tBonded residue: %s:%d -> %s" % (sel.resname[0],
                                                         sel.resid[0],
                                                         list(resnames.values())[0]))
                conect.add(other)

            # Do the renaming
            self._apply_naming_dictionary(resnames=resnames,
                                          atomnames=atomnames)

        atomsel('all').user = 1.0
        print("\n", flush=True)
        return conect

    #==========================================================================

    def _psf_to_prmtop(self):
        """
        Runs the chamber command of ParmEd to produce AMBER format input files.

        Returns:
          True if successful
        """

        # Ask the user for additional parameter files
        if self.prompt_params:
            print("\nCurrently using the following parameter files:", flush=True)
            for prm in self.parameters:
                print("  - %s" % os.path.split(prm)[1])
            print("Enter the path to the filename(s) from the current working "
                  "directory, separated by a comma, of any additional prm or str files "
                  "you wish to use.\n", flush=True)
            try:
                input = raw_input
            except NameError:
                pass
            inp = input('> ')
            if inp:
                self.parameters.extend(inp.split(','))
        else:
            print("Using the following parameter files:")
            for prm in self.parameters:
                print("  - %s" % prm.split("/")[-1])
            print("\n")

        # Begin assembling chamber input string
        args = ["-psf", self.outprefix + ".psf",
                "-crd", self.outprefix + ".pdb"]

        # Add topology and parameter arguments
        for inp in set(self.topologies + self.parameters):
            args += ["-toppar", inp]

        # Add box information since it is not in the pdb
        box = molecule.get_periodic(molid=self.molid)
        args += ["-box", " %f,%f,%f" % (box['a'], box['b'], box['c'])]
        args += ["nosettle"]

        print("Running chamber. This may take a while...", flush=True)
        parm = AmberParm()
        with warnings.catch_warnings(record=True) as w:
            action = chamber(parm, *args)
            action.execute()
            w = [i for i in w if issubclass(i.category, ParameterWarning)]

        # Do hydrogen mass repartitioning if requested
        if self.hmr:
            print("Repartitioning hydrogen masses...")
            parm = action.parm
            action = HMassRepartition(parm, "dowater")
            print(action)
            action.execute()

        print("\tRan chamber")
        write = parmout(action.parm,
                        self.outprefix + ".prmtop",
                        self.outprefix + ".inpcrd")
        write.execute()
        print("\nWrote output prmtop and inpcrd")
        return True

    #==========================================================================

    def _write_lipids(self):
        """
        Splits lipids into modular tail, head, tail that Lipid14 specifies.
        Closes the old molecule and loads the new renumbered molecule.
        Does name matching for lipids. Writes the pdb file with TER cards
        in between each lipid.

        Returns:
            (str): File name of PDB file written

        Raises:
            ValueError if an invalid lipid is found
        """
        lipid_res = set(atomsel(self.lipid_sel).residue)
        n_lips = len(lipid_res)
        if not n_lips:
            return None

        molecule.set_top(self.molid)
        f, temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_lipids_',
                                   dir=self.tmp_dir)
        fileh = os.fdopen(f, 'w')

        # Check if it's a normal residue first in case cholesterol etc in
        # the selection
        resid = 1
        idx = 1
        while lipid_res:
            residue = lipid_res.pop()
            if len(lipid_res) % 1 == 0:
                print("Writing lipids.... %.0f%%  \r"
                      % (100.-100.*len(lipid_res)/float(n_lips)), flush=True)

            sel = atomsel('residue %s' % residue)
            headres, headnam, minusidx = self.matcher.get_lipid_head(sel)

            # If it's not a lipid head, check if it's a normal residue
            if not headres:
                resnames, atomnames = self.matcher.get_names(sel, print_warning=False)
                if not resnames:
                    fileh.close()
                    raise DabbleError("Residue %s:%s not a valid lipid" %
                                      (sel.resname[0], sel.resid[0]))
                self._apply_naming_dictionary(resnames=resnames,
                                              atomnames=atomnames)
                sel.resid = resid
                resid += 1
                continue
            else:
                # Apply the name to the heads
                self._apply_naming_dictionary(resnames=headres,
                                              atomnames=headnam)

                # Pull out the tail resnames and indices
                taildicts = self.matcher.get_lipid_tails(sel, headnam.keys())
                for (resnames, atomnames) in taildicts:
                    self._apply_naming_dictionary(resnames=resnames,
                                                  atomnames=atomnames)

                # Renumber the first tail, head, then second tail and write
                # them separately. Needs to be done this way to guarantee order.
                # An atom index that's in the minus tail is given by get_lipid_head.

                # First tail
                firstdict = [_ for _ in taildicts if minusidx in _[0].keys()]
                if len(firstdict) != 1:
                    fileh.close()
                    raise DabbleError("Error finding tails for lipid %s:%s" %
                                      (sel.resname[0], sel.resid[0]))
                firstdict = firstdict[0]

                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               firstdict[0].keys()]))
                lsel.resid = resid
                lsel.user = 0.0
                idx = self._write_residue(lsel, fileh, idx)
                taildicts.remove(firstdict)

                # Head
                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               headnam.keys()]))
                lsel.resid = resid + 1
                lsel.user = 0.0
                idx = self._write_residue(lsel, fileh, idx)

                # Second tail
                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               taildicts[0][0].keys()]))
                lsel.resid = resid + 2
                lsel.user = 0.0
                idx = self._write_residue(lsel, fileh, idx)
                resid += 3
                fileh.write("TER\n") # TER card between lipid residues

        fileh.write("END\n")
        fileh.close()
        print("\n", flush=True)
        return temp

    #==========================================================================

    def _write_residue(self, ressel, fileh, idx, hetatom=True):
        """
        Writes a residue to a pdb file

        Args:
            ressel (VMD atomsel): Atom selection to write
            fileh (file handle): The file to write to
            idx (int): Current atom index
            hetatom (bool): Whether or not this is a heteroatom
        """
        for i in ressel.index:
            a = atomsel('index %d' % i) # pylint: disable=invalid-name
            fileh.write(self.get_pdb_line(a, idx, a.resid[0], hetatom=hetatom))
            idx += 1
        return idx

    #==========================================================================

    def _write_ligands(self):
        """
        Writes any remaining user=1.0 residues each to a pdb file. Renumbers
        the residues along the way. Previously this was done with a mol2 file,
        but tleap takes charges from there and we will use a topology library
        file to get the connectivity correct.

        Returns:
            (list of 2-tuple): ilename of each ligand, unit name of ligand
        """

        idx = 1
        pdbs = []

        for residue in set(atomsel("user 1.0").residue):
            _, temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_extra',
                                       dir=self.tmp_dir)
            os.close(_)
            sel = atomsel("residue %d" % residue)
            sel.user = 0.0
            sel.resid = idx
            sel.write("pdb", temp)
            unit = self.matcher.get_unit(sel)
            pdbs.append((temp, unit))
            idx += 1
        return pdbs

    #==========================================================================

    def _write_solvent(self):
        """
        Writes ions. Renames waters according to the water model and writes
        multiple PDB files to bypass PDB file format issues.

        Returns:
            (list of str): PDB filenames that were written
        """
        # Match up ion names to topology templates
        written = []
        allions = []
        for resname in set(atomsel("numbonds 0").resname):
            allions.extend(self._rename_by_resname(resname, renumber=True))

        # Dump to a temporary PDB file
        if allions:
            _, temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_ion',
                                       dir=self.tmp_dir)
            os.close(_)
            ionsel = atomsel("residue %s" % ' '.join(str(_) for _ in allions))
            ionsel.write('pdb', temp)
            ionsel.user = 0.0
            written.append(temp)

        # Now write water PDBs, waters are already named
        written.extend(self._write_water_pdbs())

        return written

    #==========================================================================

    def _write_protein(self):
        """
        Writes the protein chain to a pdb file. Writes each fragment to a
        separate file. This is somewhat inelegant but necessary in order to
        have disulfide bonds between different fragments have sane residue
        numbering in leap.

        Returns:
            list of (str, str) Name of pdb files written, and UNIT sequence
                inside
        """
        pdbinfos = []
        psel = atomsel("user 1.0 and resname %s"
                       % " ".join(self.matcher.AMINO_ACIDS))

        # Start the protein numbering from 1, as it's lost in the prmtop anyway.
        # Make resids increase across chains/fragments too, so that bond section
        # for covalent modifications is simple.
        fragres = 1
        for i, frag in enumerate(sorted(set(psel.fragment))):
            f, temp = tempfile.mkstemp(suffix='_prot.pdb', prefix='amber_prot_',
                                       dir=self.tmp_dir)

            # Now write out all the resides to a pdb file
            resseq = []
            with os.fdopen(f, 'w') as fileh:
                idx = 1
                # Grab resids again since they may have updated
                # Check for multiple residues with the same resid (insertion codes)
                for resid in sorted(set(atomsel("fragment '%s'" % frag).resid)):
                    selstr = "fragment '%s' and resid '%d' " \
                             "and user 1.0" % (frag, resid)
                    for residue in sorted(set(atomsel(selstr).residue)):
                        sel = atomsel("residue %d" % residue)
                        sel.resid = fragres
                        sel.segname = str(i)
                        sel.user = 0.0
                        idx = self._write_residue(sel, fileh, idx, hetatom=False)
                        resseq.append(sel.resname[0])
                        fragres += 1
                fileh.write("END\n")
                fileh.close()

            # Insert line breaks every 100 AA to avoid buffer overflow in tleap
            seqstring = ""
            for _, res in enumerate(resseq):
                seqstring += " %s" % res
                if _ % 100 == 0:
                    seqstring += "\n"
            pdbinfos.append((temp, seqstring))

        return pdbinfos

    #==========================================================================

    def _run_leap(self, ligfiles, prot_pdbseqs, pdbs, conect):
        """
        Runs leap, creating a prmtop and inpcrd from the given pdb and off
        library files.

        Args:
            ligfiles (dict str -> str): UNIT name and filename of mol2 file
                for each ligand. The unit name is necessary here to add
                the right variable names in leap because it is the worst.
            prot_pdbseq (tuple str,str): PDB file containing protein fragments,
                sequence of UNITs for those fragments
            pdbs (list of str): PDB or Mol2 files to combine
            conect (set of int): Atom indices connected by an extraresidue bond

        Returns:
            (str) Prefix of file written

        Raises:
            DabbleError if tleap is not present
            DabbleError if topology type cannot be determined
        """
        # Ensure leap is actually available
        tleap = os.path.join(os.environ.get("AMBERHOME"), "bin", "tleap")
        if not (os.path.isfile(tleap) and os.access(tleap, os.X_OK)):
            raise DabbleError("$AMBERHOME/bin/tleap not found at '%s'. "
                              "Check your amber installation." % tleap)

        # Create the leap input file
        f, leapin = tempfile.mkstemp(suffix='.in', prefix='dabble_leap_',
                                     dir=self.tmp_dir)
        with os.fdopen(f, 'w') as fileh:
            for i in self.topologies + self.parameters:
                if "leaprc" in i:
                    fileh.write("source %s\n" % i)
                elif "frcmod" in i:
                    fileh.write("loadamberparams %s\n" % i)
                elif ".lib" in i:
                    fileh.write("loadoff %s\n" % i)
                elif ".off" in i:
                    continue
                else:
                    raise DabbleError("Unknown topology type: %s" % i)
            fileh.write('\n')

            # Alias for water resnames (already done in water leaprc, but
            # we do it here explicitly). This makes visualization easier since
            # the water resname is wat
            fileh.write("WAT = %s\n" % self.WATER_NAMES[self.water_model])

            # Add off files here
            for i in [_ for _ in self.topologies + self.parameters if ".off" in _]:
                fileh.write("loadoff %s\n" %i)

            pdbs = [_ for _ in pdbs if _ is not None]
            for i, pdb in enumerate(pdbs):
                if "pdb" in pdb:
                    fileh.write("p%s = loadpdb %s\n" % (i, pdb))
                elif "mol2" in pdb:
                    fileh.write("p%s = loadmol2 %s\n" % (i, pdb))
                else:
                    raise DabbleError("Unknown coordinate type: %s"
                                      % pdb)

            for i, f in enumerate(ligfiles):
                if "pdb" in f[0]:
                    fileh.write("l%s = loadpdbusingseq %s {%s}\n"
                                % (i, f[0], f[1]))
                elif "mol2" in f[0]:
                    fileh.write("l%s = loadmol2 %s\n" % (i, f[0]))
                else:
                    raise DabbleError("Unknown ligand file type: %s" % f[0])

            for i, pp in enumerate(prot_pdbseqs):
                fileh.write("pp%d = loadpdbusingseq %s { %s} \n" % (i, pp[0], pp[1]))

            # Need to combine before creating bond lines since can't create
            # bonds between UNITs
            if prot_pdbseqs:
                fileh.write("p = combine { %s }\n"
                            % ' '.join(["pp%d" % i
                                        for i in range(len(prot_pdbseqs))]))

            # Create bond lines
            while conect:
                # Pull out two atoms bound to each other
                idx = conect.pop()
                s1 = atomsel("index %d" % idx)
                other = [s for s in s1.bonds[0] if s in conect]
                if len(other) != 1:
                    raise ValueError("Problem with bonds to index %d" % idx)
                other = other[0]
                s2 = atomsel("index %d" % other)
                conect.remove(other)

                fileh.write("bond p.{0}.{1} p.{2}.{3}\n".format(s1.resid[0],
                                                                s1.name[0],
                                                                s2.resid[0],
                                                                s2.name[0]))

            if pdbs:
                fileh.write("\np = combine { p %s }\n"
                            % ' '.join(["p%d"%i for i in range(len(pdbs))]))
            if ligfiles:
                fileh.write("p = combine { p %s }\n"
                            % ' '.join(["l%d"%i for i in range(len(ligfiles))]))
            fileh.write("setbox p centers 0.0\n")
            fileh.write("saveamberparm p %s.prmtop %s.inpcrd\n"
                        % (self.outprefix, self.outprefix))
            fileh.write("quit\n")

        # Now invoke leap. If it fails, print output
        out = ""
        try:
            out = check_output([tleap, "-f", leapin]).decode("utf-8")
            out = "\n===============BEGIN TLEAP OUTPUT===============\n" + out \
                + "\n================END TLEAP OUTPUT================\n"

            if self.debug:
                print(out)
            if "not saved" in out:
                raise DabbleError("Tleap call failed")
        except:
            print(out)
            raise DabbleError("Call to tleap failed! See above output for "
                              "more information")

        # Do a quick sanity check that all the protein is present.
        mademol = molecule.load("parm7", self.outprefix+ ".prmtop",
                                "rst7", self.outprefix+ ".inpcrd")
        checkacids = " ".join(self.matcher.AMINO_ACIDS)
        if len(atomsel("resname %s" % checkacids, molid=mademol)) \
           != len(atomsel("resname %s" % checkacids, molid=self.molid)):
            print(out)
            raise DabbleError("Not all protein was present in the output "
                              "prmtop. This indicates a problem with tleap. "
                              "Check the above output, especially for covalent "
                              "ligands. Is naming consistent in all .off "
                              "files?")

        return self.outprefix

    #==========================================================================

    def _alter_leaprcs(self):
        """
        Makes a copy of all passed leaprcs. Removes lines that will cause
        leap to make assumptions about the system that could break it, such
        as pdb terminal residues.
        """

        new_topos = []
        for t in self.topologies:
            if "leaprc" not in t:
                new_topos.append(t)
                continue

            f, temp = tempfile.mkstemp(suffix=".leaprc", prefix="amber_",
                                       dir=self.tmp_dir)
            inpdb = False
            with open(t, 'r') as fn:
                with os.fdopen(f, 'w') as tf:
                    for line in fn:
                        if inpdb and line == "}\n":
                            inpdb = False
                        elif "addpdbresmap" in line.lower():
                            inpdb = True

                        if not inpdb:
                            tf.write(line)

                    tf.close()
                    fn.close()
            new_topos.append(temp)

        self.topologies = new_topos

    #========================================================================#
    #                            Static methods                              #
    #========================================================================#

    @classmethod
    def get_topologies(cls, forcefield, water_model):

        if forcefield == "amber":
            if not os.environ.get("AMBERHOME"):
                raise DabbleError("AMBERHOME must be set to use AMBER "
                                  "forcefield!")

        # Check amber version and set topologies accordingly
        ambpath = os.path.join(os.environ["AMBERHOME"], "dat", "leap", "cmd")
        topologies = [
            "leaprc.protein.ff14SB",
            "leaprc.lipid14",
            "leaprc.gaff2",
        ]
        if water_model == "tip3":
            topologies.append("leaprc.water.tip3p")
        elif water_model == "tip4e":
            topologies.append("leaprc.water.tip4pew")
        elif water_model == "spce":
            topologies.append("leaprc.water.spce")
        else:
            raise DabbleError("Water model '%s' not supported with AMBER"
                              % water_model)

        for i, top in enumerate(topologies):
            topologies[i] = os.path.abspath(os.path.join(ambpath, top))
            if not os.path.isfile(topologies[i]):
                raise DabbleError("AMBER forcefield files '%s' not found\n"
                                  "Dabble requires >= AmberTools16" % top)
        return topologies

    #========================================================================#

    @classmethod
    def get_parameters(cls, forcefield, water_model):
        # AMBER topologies and parameters use the same leaprcs
        return cls.get_topologies(forcefield, water_model)

    #========================================================================#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
