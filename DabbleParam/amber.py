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
from pkg_resources import resource_filename

import vmd
import molecule
from atomsel import atomsel

from DabbleParam import CharmmWriter
from parmed.tools import chamber, parmout, HMassRepartition, defineSolvent
from parmed.amber import AmberParm

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
        self.topologies = ''
        if self.forcefield is 'charmm':
            self.parameters = [
                resource_filename(__name__, "charmm_parameters/toppar_water_ions.str"),
                resource_filename(__name__, "charmm_parameters/par_all36_cgenff.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_prot.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_lipid.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_carb.prm"),
                resource_filename(__name__, "charmm_parameters/par_all36_na.prm"),
                resource_filename(__name__, "charmm_parameters/toppar_all36_prot_na_combined.str")
                ]
        else:
            self.parameters = [] # TODO amber default parameters

        if extra_params is not None:
            self.parameters.extend(extra_params)
            self.prompt_params = False
        else:
            self.prompt_params = True

    #==========================================================================

    def write(self, prmtop_name):
        """
        Creates a prmtop with either AMBER or CHARMM parameters.
        """
        self.prmtop_name = prmtop_name

        if self.forcefield is 'charmm':
            psfgen = CharmmWriter(molid=self.molid, 
                                  tmp_dir=self.tmp_dir,
                                  lipid_sel=self.lipid_sel,
                                  extra_topos=self.extra_topos)
            self.topologies = psfgen.write(self.prmtop_name)
            self._psf_to_charmm_amber()
        elif self.forcefield is 'amber':
            # TODO TODO TODO
            print("Currently unsupported")
            quit(1)

    #========================================================================#
    #                           Private methods                              #
    #========================================================================#

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
        for top in self.topologies:
            args += ' -top %s' % top
        for prm in self.parameters:
            if 'toppar' in prm:
                args += ' -toppar %s' % prm
            else:
                args += ' -param %s' % prm

        # Add box information since it is not in the pdb
        box = molecule.get_periodic(molid=self.molid)
        args += " -box %f,%f,%f" % (box['a'], box['b'], box['c'])

        print("\nINFO: Running chamber. This may take a while...")
        sys.stdout.flush()
        parm = AmberParm()
        action = chamber(parm, args)
        print(action)
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


