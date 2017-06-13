"""
 This module contains the AmberWriter class and associated methods,
 which outputs a prmtop/inpcrd file with CHARMM names and parameters
 for use with simulation with the AMBER molecular dynamics package.
 It does this by using the chamber functionality of ParmEd API.

 Author: Robin Betz

 Copyright (C) 2015 Robin Betz
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
import sys
import tempfile
import warnings
from subprocess import check_output
from pkg_resources import resource_filename
from vmd import molecule, atomsel

from networkx.drawing.nx_pydot import write_dot
from parmed.tools import chamber, parmout, HMassRepartition, checkValidity
from parmed.amber import AmberParm
from parmed.exceptions import ParameterWarning
from Dabble.molutils import DabbleError
from Dabble.param import CharmmWriter, AmberMatcher

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberWriter(object):
    """
    Creates AMBER-format input files, using either CHARMM or AMBER parameters.

    When using the CHARMM parameters, creates a psf/pdb file and interfaces
    with ParmEd chamber command to create AMBER format files.
    When using the AMBER parameters, creates a pdb file and runs the amber
    lipid conversion script to create leap input files.
    """

    #==========================================================================

    def __init__(self, molid, **kwargs):
        """
        Creates an AMBER Writer.

        Args:
            molid (int): VMD molecule ID of system to write
            tmp_dir (str): Directory for temporary files. Defaults to "."
            forcefield (str): charmm36mm, charmm36, or amber
            lipid_sel (str): Lipid selection string. Defaults to "lipid"
            hmr (bool): If hydrogen masses should be repartitioned. Defaults
                to False.
            extra_topos (list of str): Additional topology (.str, .off, .lib) to
                include.
            extra_params (list of str): Additional parameter sets (.str, .frcmod)
            override_defaults (bool): If set, omits default amber ff14 parameters.
            debug_verbose (bool): Prints additional output, like from tleap.

        """

        self.molid = molid
        self.prmtop_name = ""

        self.tmp_dir = kwargs.get("tmp_dir", ".")
        self.lipid_sel = kwargs.get("lipid_sel", "lipid")
        self.hmr = kwargs.get("hmr", False)
        self.extra_topos = kwargs.get("extra_topos", None)
        self.override = kwargs.get("override_defaults", False)
        self.debug_verbose = kwargs.get("debug_verbose", False)

        forcefield = kwargs.get("forcefield", "charmm36m")
        if forcefield not in ["amber", "charmm36m", "charmm", "charmm36"]:
            raise DabbleError("Unsupported forcefield: %s" % forcefield)
        self.forcefield = forcefield
        if self.forcefield == "charmm36m":
            self.parameters = [
                resource_filename(__name__, "charmm_parameters/toppar_water_ions.str"),
                resource_filename(__name__, "charmm_parameters/par_all36_cgenff.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36m_prot.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_lipid.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_carb.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_na.prm"),
                resource_filename(__name__, "charmm_parameters/toppar_all36_prot_na_combined.str")
                ]
            self.topologies = []
        elif self.forcefield in ["charmm36", "charmm"]:
            self.parameters = [
                resource_filename(__name__, "charmm_parameters/toppar_water_ions.str"),
                resource_filename(__name__, "charmm_parameters/par_all36_cgenff.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_prot.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_lipid.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_carb.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_na.prm"),
                resource_filename(__name__, "charmm_parameters/toppar_all36_prot_na_combined.str")
                ]
            self.topologies = []
        elif self.forcefield == 'amber':
            if not os.environ.get("AMBERHOME"):
                raise DabbleError("AMBERHOME must be set to use AMBER forcefield!")
            if not os.path.isfile(os.path.join(os.environ.get("AMBERHOME"),
                                               "bin", "tleap")):
                raise DabbleError("tleap is not present in $AMBERHOME/bin!")

            # Check amber version and set topologies accordingly
            self.topologies = [
                "leaprc.protein.ff14SB",
                "leaprc.lipid14",
                "leaprc.water.tip3p",
                "leaprc.gaff",
            ]
            for i, top in enumerate(self.topologies):
                self.topologies[i] = os.path.join(os.environ["AMBERHOME"],
                                                  "dat", "leap", "cmd", top)
                if not os.path.isfile(self.topologies[i]):
                    raise DabbleError("AMBER version too old! "
                                      "Dabble requires >= AmberTools16!")

            self.parameters = []
            self.matcher = None

        if self.override:
            self.topologies = []
            self.parameters = []

        if kwargs.get("extra_topos") is not None:
            self.topologies.extend(kwargs.get("extra_topos"))

        if kwargs.get("extra_params") is not None:
            self.parameters.extend(kwargs.get("extra_params"))

        self.prompt_params = False

    #==========================================================================

    def write(self, prmtop_name):
        """
        Creates a prmtop with either AMBER or CHARMM parameters.
        """
        self.prmtop_name = prmtop_name

        # Charmm forcefield
        if "charmm" in self.forcefield:
            psfgen = CharmmWriter(molid=self.molid,
                                  tmp_dir=self.tmp_dir,
                                  lipid_sel=self.lipid_sel,
                                  extra_topos=self.extra_topos,
                                  override_defaults=self.override)
            self.topologies = psfgen.write(self.prmtop_name)
            self._psf_to_charmm_amber()

        # Amber forcefield
        elif "amber" in self.forcefield:

            # Initialize the matcher
            self.matcher = AmberMatcher(self.topologies)
            print("Using the following topologies:")
            for top in self.topologies:
                print("  - %s" % top.split("/")[-1])
                top = os.path.abspath(top)

            for par in self.parameters:
                par = os.path.abspath(par)

            # Assign atom types
            print("Assigning AMBER atom types...")
            self._split_caps()
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
                write = parmout(action.parm, "%s.prmtop" % self.prmtop_name)
                                                                      #self.prmtop_name))
                write.execute()
                parm = write.parm

            # Check validity of output prmtop using parmed
            parm = AmberParm(prm_name=self.prmtop_name+".prmtop",
                             xyz=self.prmtop_name+".inpcrd")
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
            (set of tuples (int,int)): Residue #s of disulfide bonded residues

        Raises:
            ValueError if a residue definition could not be found
        """

        nonlips = set(atomsel("not (%s)" % self.lipid_sel,
                              molid=self.molid).get("residue"))
        n_res = len(nonlips)
        conect = set() # Atom indices bound to noncanonical residues
        while nonlips:
            if len(nonlips) % 500 == 0:
                sys.stdout.write("Renaming residues.... %.0f%%  \r"
                                 % (100.-100*len(nonlips)/float(n_res)))
                sys.stdout.flush()

            residue = nonlips.pop()
            sel = atomsel("residue %s" % residue)
            resnames, atomnames = self.matcher.get_names(sel, print_warning=False)

            # Check if it's a linkage to another amino acid
            if not resnames:
                resnames, atomnames, other = self.matcher.get_linkage(sel, self.molid)
                if not resnames:
                    rgraph = self.matcher.parse_vmd_graph(sel)[0]
                    write_dot(rgraph, "rgraph.dot")
                    raise DabbleError("ERROR: Could not find a residue definition "
                                      "for %s:%s" % (sel.get("resname")[0],
                                                     sel.get("resid")[0]))

                print("\tBonded residue: %s:%d -> %s" % (sel.get("resname")[0],
                                                         sel.get("resid")[0],
                                                         list(resnames.values())[0]))
                conect.add(other)

            # Do the renaming
            self._apply_naming_dictionary(resnames, atomnames)

        atomsel('all').set('user', 1.0)
        sys.stdout.write("\n")
        return conect

    #==========================================================================

    def _psf_to_charmm_amber(self):
        """
        Runs the chamber command of ParmEd to produce AMBER format input files.

        Returns:
          True if successful
        """

        # Ask the user for additional parameter files
        if self.prompt_params:
            sys.stdout.flush()
            print("\nCurrently using the following parameter files:")
            for prm in self.parameters:
                print("  - %s" % prm.split("/")[-1])
            print("Enter the path to the filename(s) from the current working "
                  "directory, separated by a comma, of any additional prm or str files "
                  "you wish to use.\n")
            sys.stdout.flush()
            try: input = raw_input
            except NameError: pass
            inp = input('> ')
            if inp:
                self.parameters.extend(inp.split(','))
        else:
            print("Using the following parameter files:")
            for prm in self.parameters:
                print("  - %s" % prm.split("/")[-1])
            print("\n")

        # Begin assembling chamber input string
        args = "-crd %s.pdb -psf %s.psf" % (self.prmtop_name, self.prmtop_name)

        # Add topology and parameter arguments
        for inp in set(self.topologies + self.parameters):
            args += ' -toppar %s' % inp

        # Add box information since it is not in the pdb
        box = molecule.get_periodic(molid=self.molid)
        args += " -box %f,%f,%f" % (box['a'], box['b'], box['c'])
        args += " nosettle"

        print("Running chamber. This may take a while...")
        sys.stdout.flush()
        parm = AmberParm()
        with warnings.catch_warnings(record=True) as w:
            action = chamber(parm, args)
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
        write = parmout(action.parm, "%s.prmtop %s.inpcrd"
                        %(self.prmtop_name, self.prmtop_name))
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
        lipid_res = set(atomsel(self.lipid_sel).get('residue'))
        n_lips = len(lipid_res)
        if not n_lips:
            return None

        molecule.set_top(self.molid)
        temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_lipids_',
                                dir=self.tmp_dir)[1]
        fileh = open(temp, 'w')

        # Check if it's a normal residue first in case cholesterol etc in
        # the selection
        resid = 1
        idx = 1
        while lipid_res:
            residue = lipid_res.pop()
            if len(lipid_res) % 1 == 0:
                sys.stdout.write("Writing lipids.... %.0f%%  \r"
                                 % (100.-100.*len(lipid_res)/float(n_lips)))
                sys.stdout.flush()

            sel = atomsel('residue %s' % residue)
            headres, headnam, minusidx = self.matcher.get_lipid_head(sel)

            # If it's not a lipid head, check if it's a normal residue
            if not headres:
                resnames, atomnames = self.matcher.get_names(sel, print_warning=False)
                if not resnames:
                    raise DabbleError("Residue %s:%s not a valid lipid" %
                                      (sel.get('resname')[0], sel.get('resid')[0]))
                self._apply_naming_dictionary(resnames, atomnames)
                sel.set('resid', resid)
                resid += 1
                continue
            else:
                # Apply the name to the heads
                self._apply_naming_dictionary(headres, headnam)

                # Pull out the tail resnames and indices
                taildicts = self.matcher.get_lipid_tails(sel, headnam.keys())
                for (resnames, atomnames) in taildicts:
                    self._apply_naming_dictionary(resnames, atomnames)

                # Renumber the first tail, head, then second tail and write
                # them separately. Needs to be done this way to guarantee order.
                # An atom index that's in the minus tail is given by get_lipid_head.

                # First tail
                firstdict = [_ for _ in taildicts if minusidx in _[0].keys()]
                if len(firstdict) != 1:
                    raise DabbleError("Error finding tails for lipid %s:%s" %
                                     (sel.get('resname')[0], sel.get('resid')[0]))
                firstdict = firstdict[0]

                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               firstdict[0].keys()]))
                lsel.set('resid', resid)
                lsel.set('user', 0.0)
                idx = self._write_residue(lsel, fileh, idx)
                taildicts.remove(firstdict)

                # Head
                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               headnam.keys()]))
                lsel.set('resid', resid+1)
                lsel.set('user', 0.0)
                idx = self._write_residue(lsel, fileh, idx)

                # Second tail
                lsel = atomsel('index %s' % ' '.join([str(x) for x in \
                               taildicts[0][0].keys()]))
                lsel.set('resid', resid+2)
                lsel.set('user', 0.0)
                idx = self._write_residue(lsel, fileh, idx)
                resid += 3
                fileh.write("TER\n") # TER card between lipid residues

        fileh.write("END\n")
        fileh.close()
        sys.stdout.write("\n")
        return temp

    #==========================================================================

    def _write_residue(self, ressel, fileh, idx, hetatm=False):
        """
        Writes a residue to a pdb file

        Args:
            ressel (VMD atomsel): Atom selection to write
            fileh (file handle): The file to write to
            idx (int): Current atom index
            hetatm (bool): Whether or not this is a heteroatom
        """
        if hetatm:
            atm = "HETATM"
        else:
            atm = "ATOM"
        for i in ressel.get('index'):
            a = atomsel('index %d' % i) # pylint: disable=invalid-name
            assert len(a) == 1
            entry = ('%-6s%5d %-5s%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f'
                     '     %-4s%2s\n' % (atm, idx, a.get('name')[0],
                                         a.get('resname')[0],
                                         a.get('chain')[0],
                                         a.get('resid')[0],
                                         a.get('x')[0],
                                         a.get('y')[0],
                                         a.get('z')[0],
                                         0.0, 0.0, a.get('segname')[0],
                                         a.get('element')[0]))
            fileh.write(entry)
            idx += 1
        return idx

    #==========================================================================

    def _apply_naming_dictionary(self, resnames, atomnames):
        """
        Applies the names from a matcher.
        """
        for idx, name in atomnames.items():
            atom = atomsel('index %s' % idx)
            if atom.get('name')[0] != name:
                atom.set('name', name)
            atom.set('resname', resnames[idx])

    #==========================================================================

    def _split_caps(self):
        """
        Pulls out ACE and NMA caps, renumbers residues, and loads that
        renumbered molecule. Closes the old molecule and sets this as
        the top one.

        Returns:
            (int): Molid of new molecule
        """
        # Put our molecule on top and grab selection
        molecule.set_top(self.molid)

        # Do one fragment at a time to handle duplicate resids across chains
        print("Checking if capping groups need to be renumbered")
        for frag in set(atomsel('all').get('fragment')):
            for resid in sorted(set(atomsel('fragment %d' % frag).get('resid'))):
                # Handle bug where capping groups in same residue as the
                # neighboring amino acid Maestro writes it this way for some
                # reason but it causes problems down the line when psfgen doesn't
                # understand the weird combined residue
                rid = atomsel("fragment '%d' and resid '%d'"
                              % (frag, resid)).get('residue')[0]
                names = set(atomsel('residue %d'% rid).get('resname'))
                assert len(names) < 3, ("More than 2 residues with same number... "
                                        "currently unhandled. Report a bug")

                if len(names) > 1:
                    if 'ACE' in names and 'NMA' in names:
                        print("ERROR: Both ACE and NMA were given the same resid"
                              "Check your input structure")
                        quit(1)

                    if 'ACE' in names:
                        # Set ACE residue number as one less
                        resid = atomsel('residue %d and not resname ACE' % rid).get('resid')[0]
                        if len(atomsel("fragment '%s' and resid '%d'" % (frag, resid-1))):
                            raise ValueError('ACE resid collision number %d' % resid-1)
                        atomsel('residue %d and resname ACE'
                                % rid).set('resid', resid-1)
                        print("\tACE %d -> %d" % (resid, resid-1))

                    elif 'NMA' in names:
                        # Set NMA residue number as one more
                        resid = atomsel('residue %d and not resname NMA' % rid).get('resid')[0]
                        if len(atomsel("fragment '%s' and resid '%d'" % (frag, resid+1))):
                            raise ValueError('NMA resid collision number %d' % resid+1)

                        atomsel('residue %d and resname NMA'
                                % rid).set('resid', resid+1)
                        print("\tNMA %d -> %d" % (resid, resid+1))

        # Have to save and reload so residues are parsed correctly by VMD
        temp = tempfile.mkstemp(suffix='.mae', prefix='mae_renum_',
                                dir=self.tmp_dir)[1]
        atomsel('all').write('mae', temp)
        molecule.delete(self.molid)
        self.molid = molecule.load('mae', temp)
        molecule.set_top(self.molid)

        return self.molid

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

        for residue in set(atomsel("user 1.0").get("residue")):
            temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_extra',
                                    dir=self.tmp_dir)[1]
            sel = atomsel("residue %d" % residue)
            sel.set('user', 0.0)
            sel.set('resid', idx)
            sel.write("pdb", temp)
            unit = self.matcher.get_unit(sel)
            pdbs.append((temp, unit))
            idx += 1
        return pdbs

    #==========================================================================

    def _write_solvent(self):
        """
        Renumbers the waters into different chains to bypass pdb file format
        issues with more than 10000 waters in a chain. Also writes ions.

        Returns:
            (list of str): PDB filenames that were written
        """

        written = []
        # First write all the ions using just vmd
        temp = tempfile.mkstemp(suffix='.pdb', prefix='amber_ion',
                                dir=self.tmp_dir)[1]
        ionnames = [_ for _ in self.matcher.known_res.keys() if '+' in _ or '-' in _]
        ions = atomsel("resname NA CL %s and user 1.0" % ' '.join("'%s'" % _ for _ in ionnames))
        if len(ions):
            ions.set('resid', range(1, len(ions)+1))
            ions.write('pdb', temp)
            ions.set('user', 0.0)
            written.append(temp)

        # Select all the unwritten waters, using the user field to track which
        # ones have been written.
        allw = atomsel('water and user 1.0')

        # Find the problem waters with unordered indices
        problems = []
        for r in set(allw.get('residue')):
            widx = atomsel('residue %s' % r).get("index")
            if max(widx) - min(widx) != 2:
                problems.append(r)
                atomsel('residue %s' % r).set("user", 0.0) # get it out of allw

        allw.update()
        num_written = int(len(allw)/(9999*3))+1

        # Pull out and write 10k waters at a time if we have normal waters
        if allw:
            print("Going to write %d files for %d water atoms"
                  % (num_written, len(allw)))
            for i in range(num_written):
                temp = tempfile.mkstemp(suffix='_%d.pdb' % i, prefix='amber_wat_',
                                        dir=self.tmp_dir)[1]
                written.append(temp)
                residues = list(set(allw.get('residue')))[:9999]

                batch = atomsel('residue %s' % ' '.join([str(x) for x in residues]))
                try:
                    batch.set('resid', [k for k in range(1, int(len(batch)/3)+1)
                                        for _ in range(3)])
                except ValueError:
                    print("\nERROR! You have some waters missing hydrogens!\n"
                          "Found %d water residues, but %d water atoms. Check "
                          " your crystallographic waters in the input structure."
                          % (len(residues), len(batch)))
                    quit(1)
                batch.set('user', 0.0)
                batch.write('pdb', temp)
                allw.update()

        # Now write the problem waters
        if problems:
            written.append(self._write_unorderedindex_waters(problems))

        return written

    #==========================================================================

    def _write_unorderedindex_waters(self, residues):
        """
        Renumbers and sorts the specified waters manually. This is much less
        efficient but is necessary in cases where atoms within a water molecule
        are not sequential in index, preventing quick renaming with VMD.
        Identify problem waters, then call this on them. It'll write its own
        psf_wat_* file with just those waters, minimizing inefficiency.

        Args:
            residues (list of int): Problem water molecules
        Returns:
            (str) Name of the pdb file written
        """
        temp = tempfile.mkstemp(suffix='_indexed.pdb', prefix='amber_wat_',
                                dir=self.tmp_dir)[1]
        fileh = open(temp, 'w')

        idx = 1
        for residx, residue in enumerate(residues):
            res = atomsel('residue %d' % residue)
            res.set('user', 0.0)

            for i in res.get('index'):
                a = atomsel('index %d' % i) # pylint: disable=invalid-name
                entry = ('%-6s%5d %-5s%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f'
                         '     %-4s%2s\n' % ('ATOM', idx, a.get('name')[0],
                                             a.get('resname')[0],
                                             a.get('chain')[0],
                                             residx+1,
                                             a.get('x')[0],
                                             a.get('y')[0],
                                             a.get('z')[0],
                                             0.0, 0.0, a.get('segname')[0],
                                             a.get('element')[0]))
                idx += 1
                fileh.write(entry)

        fileh.write('END\n')
        fileh.close()
        return temp

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
        psel = atomsel("user 1.0 and resname %s" % " ".join(self.matcher._acids))

        # Start the protein numbering from 1, as it's lost in the prmtop anyway.
        # Make resids increase across chains/fragments too, so that bond section
        # for covalent modifications is simple.
        fragres = 1
        for i, frag in enumerate(sorted(set(psel.get("fragment")))):
            temp = tempfile.mkstemp(suffix='_prot.pdb', prefix='amber_prot_',
                                    dir=self.tmp_dir)[1]

            for res in sorted(set(atomsel("fragment '%s'" % frag).get("resid"))):
                ressel = atomsel("fragment '%s' and resid '%d' and user 1.0"
                                 % (frag, res))
                ressel.set("resid", fragres)
                ressel.set("user", 2.0) # Can't select only be resid as it changes
                fragres += 1

            atomsel("fragment '%s' and user 2.0" % frag).set("user", 1.0)

            # Now write out all the resides to a pdb file
            resseq = []
            with open(temp, 'w') as fileh:
                idx = 1
                # Grab resids again since they may have updated
                for resid in sorted(set(atomsel("fragment '%s'" % frag).get("resid"))):
                    sel = atomsel("fragment '%s' and resid '%d' and "
                                  "user 1.0" % (frag, resid))
                    idx = self._write_residue(sel, fileh, idx, hetatm=False)
                    sel.set('segname', str(i))
                    sel.set('user', 0.0)
                    resseq.append(sel.get("resname")[0])
                fileh.write("END\n")
                fileh.close()

            # Insert line breaks every 100 AA to avoid buffer overflow in tleap
            seqstring = ""
            for i, res in enumerate(resseq):
                seqstring += " %s" % res
                if i % 100 == 0:
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
            ValueError if AMBERHOME is unset
            ValueError if topology type cannot be determined
        """
        # Ensure leap is actually available
        if not os.environ.get("AMBERHOME"):
            raise DabbleError("AMBERHOME must be set to use leap!")

        # Create the leap input file
        leapin = tempfile.mkstemp(suffix='.in', prefix='dabble_leap_',
                                  dir=self.tmp_dir)[1]
        with open(leapin, 'w') as fileh:
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
            fileh.write("p = combine { %s }\n"
                         % ' '.join(["pp%d" % i for i in range(len(prot_pdbseqs))]))

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

                fileh.write("bond p.{0}.{1} p.{2}.{3}\n".format(
                    s1.get('resid')[0],
                    s1.get('name')[0],
                    s2.get('resid')[0],
                    s2.get('name')[0]))

            if len(pdbs):
                fileh.write("\np = combine { p %s }\n"
                            % ' '.join(["p%d"%i for i in range(len(pdbs))]))
            if len(ligfiles):
                fileh.write("p = combine { p %s }\n"
                            % ' '.join(["l%d"%i for i in range(len(ligfiles))]))
            fileh.write("setbox p centers 0.0\n")
            fileh.write("saveamberparm p %s.prmtop %s.inpcrd\n"
                        % (self.prmtop_name, self.prmtop_name))
            fileh.write("quit\n")
            fileh.close()

        # Now invoke leap. If it fails, print output
        out = ""
        try:
            out = check_output([os.path.join(os.environ.get("AMBERHOME"),
                                             "bin", "tleap"),
                                "-f", leapin]).decode("utf-8")
            out = "%s%s%s" % (
                "\n================BEGIN TLEAP OUTPUT================\n",
                out,
                "\n=================END TLEAP OUTPUT=================\n")

            if self.debug_verbose:
                print(out)
            if "not saved" in out:
                raise DabbleError("Tleap call failed")
        except:
            print(out)
            raise DabbleError("Call to tleap failed! See above output for errors")
            quit(1)

        # Do a quick sanity check that all the protein is present.
        mademol = molecule.load("parm7", "%s.prmtop" % self.prmtop_name,
                                "rst7", "%s.inpcrd" % self.prmtop_name)
        if len(atomsel("resname %s" % " ".join(self.matcher._acids), mademol)) \
                != len(atomsel("resname %s" % " ".join(self.matcher._acids), self.molid)):
           print(out)
           raise DabbleError("Not all protein was present in the output prmtop."
                             " This indicates a problem with tleap. Check the "
                             "above output, especially for covalent ligands. "
                             "Is naming consistent in all .off files?")

        return self.prmtop_name

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

            temp = tempfile.mkstemp(suffix=".leaprc", prefix="amber_",
                                    dir=self.tmp_dir)[1]
            inpdb = False
            with open(t, 'r') as fn:
                with open(temp, 'w') as tf:
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


