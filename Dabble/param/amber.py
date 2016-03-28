"""
This module contains the AmberWriter class and associated methods,
which outputs a prmtop/inpcrd file with CHARMM names and parameters
for use with simulation with the AMBER molecular dynamics package.
It does this by using the chamber functionality of ParmEd API.

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
import os
import sys
import tempfile
from pkg_resources import resource_filename

import vmd
import molecule
from atomsel import atomsel

from parmed.tools import chamber, parmout, HMassRepartition
from parmed.amber import AmberParm
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

    def __init__(self, molid, tmp_dir,
                 forcefield='charmm', lipid_sel="lipid",
                 hmr=False, extra_topos=None, extra_params=None):
        self.lipid_sel = lipid_sel
        self.molid = molid
        self.tmp_dir = tmp_dir
        self.hmr = hmr
        self.extra_topos = extra_topos
        if forcefield not in ['amber', 'charmm']:
            raise ValueError("Unsupported forcefield: %s" % forcefield)
        self.forcefield = forcefield
        if self.forcefield == 'charmm':
            self.parameters = [
                resource_filename(__name__, "charmm_parameters/toppar_water_ions.str"),
                resource_filename(__name__, "charmm_parameters/par_all36_cgenff.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_prot.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_lipid.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_carb.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_na.prm"),
                resource_filename(__name__, "charmm_parameters/toppar_all36_prot_na_combined.str")
                ]
            self.topologies = ''
        elif self.forcefield == 'amber':
            if not os.environ.get("AMBERHOME"):
                raise ValueError("AMBERHOME must be set to use AMBER forcefield!")

            self.topologies = [
                os.path.join(os.environ["AMBERHOME"],"dat","leap","cmd","leaprc.ff14SB"),
                os.path.join(os.environ["AMBERHOME"],"dat","leap","cmd","leaprc.lipid14"),
            #    os.path.join(os.environ["AMBERHOME"],"dat","leap","cmd","leaprc.gaff")
               ]
            self.parameters = self.topologies

        if extra_params is not None:
            self.parameters.extend(extra_params)

        self.prompt_params = False

    #==========================================================================

    def write(self, prmtop_name):
        """
        Creates a prmtop with either AMBER or CHARMM parameters.
        """
        self.prmtop_name = prmtop_name

        # Charmm forcefield
        if self.forcefield is 'charmm':
            psfgen = CharmmWriter(molid=self.molid, 
                                  tmp_dir=self.tmp_dir,
                                  lipid_sel=self.lipid_sel,
                                  extra_topos=self.extra_topos)
            self.topologies = psfgen.write(self.prmtop_name)
            self._psf_to_charmm_amber()

        # Amber forcefield
        elif self.forcefield is 'amber':
            ## Save and reload so residue looping is correct
            self._split_caps()
            self._rename_atoms_amber()

    #========================================================================#
    #                           Private methods                              #
    #========================================================================#
    
    def _rename_atoms_amber(self):
        """
        Matches up atom names with those in the provided topologies and
        sets the atom and residue names correctly in the built molecule.

        Raises:
            ValueError if a residue definition could not be found
        """

        matcher = AmberMatcher(self.topologies)

        for residue in set(atomsel("all", molid=self.molid).get("residue")):
            sel = atomsel("residue %s" % residue)
            resnames, atomnames = matcher.get_names(sel, print_warning=False)
            
            # Check if it's disulfide bond
            if not resnames:
                resnames, atomnames, conect = matcher.get_disulfide(sel, self.molid)
                if not resnames:
                    import networkx as nx
                    rgraph = matcher.parse_vmd_graph(sel)[0]
                    nx.write_dot(rgraph, "rgraph.dot")
                    raise ValueError("ERROR: Could not find a residue definition "
                                     "for %s:%s" % (sel.get("resname")[0],
                                                    sel.get("resid")[0]))

            # Do the renaming
            for idx, name in atomnames.iteritems():
                atom = atomsel('index %s' % idx)
                if atom.get('name')[0] != name:
                    #print("Renaming %s:%s: %s -> %s:%s" % (atom.get('resname')[0], residue,
                    #                                          atom.get('name')[0],
                    #                                          resnames[idx], name))
                    atom.set('name', name)
                atom.set('resname', resnames[idx])

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
            inp = raw_input('> ')
            if inp:
                self.parameters.extend(inp.split(','))
        else:
            print("Using the following parameter files:")
            for prm in self.parameters:
                print("  - %s" % prm.split("/")[-1])

        # Begin assembling chamber input string
        args = "-crd %s.pdb -psf %s.psf" % (self.prmtop_name, self.prmtop_name)

        # Add topology and parameter arguments
        for inp in self.topologies + self.parameters:
            args += ' -toppar %s' % inp

        # Add box information since it is not in the pdb
        box = molecule.get_periodic(molid=self.molid)
        args += " -box %f,%f,%f" % (box['a'], box['b'], box['c'])

        print("\nINFO: Running chamber. This may take a while...")
        sys.stdout.flush()
        parm = AmberParm()
        action = chamber(parm, args)
        action.execute()

        # Do hydrogen mass repartitioning if requested
        if self.hmr:
            print("\nINFO: Repartitioning hydrogen masses...")
            parm = action.parm
            action = HMassRepartition(parm, "dowater")
            print(action)
            action.execute()

        print("\nINFO: Ran chamber")
        write = parmout(action.parm, "%s.prmtop %s.inpcrd"
                        %(self.prmtop_name, self.prmtop_name))
        write.execute()
        print("\nINFO: Wrote output prmtop and inpcrd")
        return True


    #==========================================================================

    def _write_amber_pdb(self):
        """
        Writes a pdb file of the selected molecule that's compatible with leap

        Args:
          output_filename (str): Name of file to save, including .pdb extension
          molid (int): VMD index of the molecule to save

        Returns:
          True if successful
        """
#    import charmmlipid2amber.py

#    oh hi write pdb here pls
#    TODO do atom names need to be changed? - no not from charmm names
#    Looks like I can save pdb/psf and then call charmmlipid2amber.py on  pdb
# oh also call leap kthx

    #==========================================================================

    def _split_caps(self):
        """
        Pulls out ACE and NMA caps, renumbers residues, and loads that 
        renumbered molecule. Closes the old molecule and sets this as
        the top one.
        
        Returns:
            (int): Molid of loaded fragment
        """
        # Put our molecule on top and grab selection
        molecule.set_top(self.molid)

        # Do one fragment at a time to handle duplicate resids across chains
        for frag in set(atomsel('all').get('fragment')):
            for resid in sorted(set(atomsel('fragment %d' % frag).get('resid'))):
                # Handle bug where capping groups in same residue as the
                # neighboring amino acid Maestro writes it this way for some
                # reason but it causes problems down the line when psfgen doesn't
                # understand the weird combined residue
                rid = atomsel('fragment %d and resid %d' % (frag, resid)).get('residue')[0]
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
                        if len(atomsel('fragment %s and resid %d' % (frag, resid-1))):
                            raise ValueError('ACE resid collision number %d' % resid-1)
                        atomsel('residue %d and resname ACE'
                                % rid).set('resid', resid-1)
                        print("INFO: ACE %d -> %d" % (resid, resid-1))

                    elif 'NMA' in names:
                        # Set NMA residue number as one more
                        resid = atomsel('residue %d and not resname NMA' % rid).get('resid')[0]
                        if len(atomsel('fragment %s and resid %d' % (frag, resid+1))):
                            raise ValueError('NMA resid collision number %d' % resid+1)

                        atomsel('residue %d and resname NMA'
                                % rid).set('resid', resid+1)
                        print("INFO: NMA %d -> %d" % (resid, resid+1))

        # Have to save and reload so residues are parsed correctly by VMD
        temp = tempfile.mkstemp(suffix='.mae', prefix='mae_renum_',
                                dir=self.tmp_dir)[1]
        atomsel('all').write('mae', temp)
        molecule.delete(self.molid)
        self.molid = molecule.load('mae', temp)
        molecule.set_top(self.molid)

        return self.molid

    #==========================================================================

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


