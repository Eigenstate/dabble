# System builder
#
# Author: Robin Betz
#
# Copyright (C) 2015 Robin Betz
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
from pkg_resources import resource_filename
import numpy as np
import math
import random
import os
import tempfile

from dabble import fileutils, molutils, DabbleError
from vmd import atomsel, molecule, trans

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants
_MEMBRANE_HYDROPHOBIC_THICKNESS = 30.0
_MEMBRANE_FULL_THICKNESS = 50.0

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DabbleBuilder(object):
    """
    Builds protein systems, optionally in a membrane.

    Tiles box appropriately, inserts protein into box,
    removes conflicts, trims excess solvent and adds ions.

    Attributes:
      molids (dict str->int): Molecule IDs comprising system components
      size (array len 3 of floats): Size of the x, y, and z dimensions of the system
      solute_sel (str): VMD atom selection string for original solute
      opts (dictionary): All options passed to the system builder
      tmp_dir (str): Directory in which to save temporary files
      water_only (bool): If the solvent is just a water box
    """

    #==========================================================================

    def __init__(self, **kwargs):
        random.seed(2015) # Be deterministic
        self.opts = kwargs
        self.molids = {}
        self.size = [0., 0., 0.]
        self.solute_sel = ""
        self.water_only = False
        self._zmax = self._zmin = 0.
        if self.opts.get('tmp_dir'):
            self.tmp_dir = self.opts.get('tmp_dir')
            if not os.path.exists(self.tmp_dir):
                os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='dabble', dir=os.getcwd())

        # Check for default lipid membrane
        # Set some default options
        if 'forcefield' not in self.opts:
            self.opts['forcefield'] = "charmm"
        if 'lipid_sel' not in self.opts:
            self.opts['lipid_sel'] = "lipid or resname POPS POPG"
        if 'cation' not in self.opts:
            self.opts['cation'] = 'Na'
        if 'anion' not in self.opts:
            self.opts['anion'] = 'Cl'
        if 'salt_conc' not in self.opts:
            self.opts['salt_conc'] = 0.150
        if 'wat_buffer' not in self.opts:
            self.opts['wat_buffer'] = 20.0
        if 'xy_buf' not in self.opts:
            self.opts['xy_buf'] = 17.5
        if 'lipid_dist' not in self.opts:
            self.opts['lipid_dist'] = 1.75
        if 'membrane_system' not in self.opts:
            self.opts['membrane_system'] = \
                resource_filename(__name__, os.path.join("lipid_membranes",
                                                         "popc.mae"))

        # Process input arguments
        if self.opts.get('membrane_system') == 'water':
            self.opts['membrane_system'] = resource_filename(__name__, \
                    "lipid_membranes/tip3pbox.mae")


        # Check the file output format is supported
        self.out_fmt = fileutils.check_out_type(
            self.opts.get('output_filename'),
            outformat=self.opts.get('format'),
            forcefield=self.opts.get('forcefield'),
            hmr=self.opts.get('hmassrepartition')
        )

    #==========================================================================

    def _build(self):
        """
        Builds the system

        Returns:
         (int) VMD molecule id of built system
        """
        # Load the membrane system, check if it's water only
        self.add_molecule(self.opts['membrane_system'], 'membrane')
        if not atomsel(self.opts['lipid_sel'], molid=self.molids['membrane']):
            self.water_only = True
            print("No lipid detected. Proceeding with pure liquid solvent")
        dims = molutils.get_system_dimensions(molid=self.molids['membrane'])
        print("Solvent patch dimensions are %.2f x %.2f x %.2f" % dims)

        # Load the solute (protein or ligand)
        print("Loading and orienting the solute...")
        self.add_molecule(self.opts['solute_filename'], 'solute')
        self._set_solute_sel(self.molids['solute'])

        # Orient the solute in the X,Y, and optionally Z directions
        self.molids['solute'] = self._orient_solute(self.molids['solute'])

        # Compute dimensions of the input system
        print("Computing the size of the input periodic cell...")
        dx_sol, dy_sol, dx_tm, dy_tm, dz_full = \
                self.get_cell_size(mem_buf=self.opts.get('xy_buf'),
                                   wat_buf=self.opts.get('wat_buffer'),
                                   molid=self.molids['solute'])
        if self.water_only:
            print("\tSolute x diameter is %.2f\n"
                  "\tSolute y diameter is %.2f\n"
                  "\tSolute z diameter is %.2f"
                  % (dx_sol, dy_sol, dz_full))
        else:
            print("\tSolute x diameter is %.2f (%.2f in TM region)\n"
                  "\tSolute y diameter is %.2f (%.2f in TM region)\n"
                  "\tSolute + membrane Z diameter is %.2f"
                  % (dx_sol, dx_tm, dy_sol, dy_tm, dz_full))

        # Compute dimensions of the final system
        print("Final system will be %.2f x %.2f x %.2f"
              % (self.size[0], self.size[1], self.size[2]))
        print("\tX,Y solvent buffer: (%4.1f, %4.1f)" % (self.size[0] - dx_sol,
                                                        self.size[1] - dy_sol))
        if not self.water_only:
            print("\tX,Y transmembrane buffer: (%4.1f, %4.1f)"
                  % (self.size[0] - dx_tm, self.size[1] - dy_tm))
        print("\tZ solvent buffer: %4.1f" % (self.size[2]-dz_full))

        # Tile the membrane/solvent
        print("Tiling solvent...")
        self.molids['tiled_membrane'], times = \
                tile_membrane_patch(self.molids['membrane'],
                                    self.size,
                                    tmp_dir=self.tmp_dir,
                                    allow_z_tile=self.water_only)
        print("\tSolvent tiled %d x %d x %d times" % times)

        # Only delete if a new molecule was created (if tiling occured)
        if self.molids['tiled_membrane'] != self.molids['membrane']:
            self.remove_molecule('membrane')

        print("Centering solvent...")
        self.molids['tiled_membrane'] = \
                molutils.center_system(molid=self.molids['tiled_membrane'],
                                       tmp_dir=self.tmp_dir, center_z=True)

        # Combine tiled membrane with solute
        print("Combining solute and tiled solvent patch...")
        self.molids['inserted'] = \
                molutils.combine_molecules(input_ids=[self.molids['solute'],
                                                      self.molids['tiled_membrane']],
                                           tmp_dir=self.tmp_dir)
        self.remove_molecule('tiled_membrane')
        self.remove_molecule('solute')

        # Add more waters if necessary
        self.molids['combined'] = self._add_water(self.molids['inserted'])
        self.remove_molecule('inserted')

        # Remove atoms outside the final system cell, accounting for boundary
        self._set_cell_to_square_prism(self.molids['combined'])
        box_wat = self._trim_water(self.molids['combined'])

        # Remove lipids outside the box if there are lipids
        if not self.water_only:
            print("\nInitial membrane composition:\n%s" %
                  molutils.print_lipid_composition(self.opts.get('lipid_sel'),
                                                   molid=self.molids['combined']))

        # Remove atoms that clash with the solvent
        clashes = self._remove_overlapping_residues(self.opts.get('lipid_sel'),
                                                    self.molids['combined'],
                                                    self.opts.get('lipid_friendly_sel'))

        print("Removed %d extra atoms\n"
              "\t%d out of the box\n"
              "\t%d too close to the solute" % ((box_wat+clashes), box_wat,
                                                clashes))

        # Remove extra lipids and print info about membrane
        if not self.water_only:
            self._remove_clashing_lipids(self.molids['combined'],
                                         self.opts.get('lipid_sel'),
                                         self.opts.get('lipid_friendly_sel'))
            print("\nFinal membrane composition:\n%s" %
                  molutils.print_lipid_composition(self.opts.get('lipid_sel'),
                                                   self.molids['combined']))

        # Calculate charge
        print("\nSolute net charge: %+d\n"
              "System net charge: %+d"
              % (molutils.get_net_charge(self.solute_sel, self.molids['combined']),
                 molutils.get_system_net_charge(self.molids['combined'])))

        # Add ions as necessary
        self.convert_ions(self.opts.get('salt_conc'),
                          molid=self.molids['combined'])

        # System is now built
        return self.molids['combined']

    #==========================================================================

    def write(self):
        """
        Writes the final built system.

        Args:
          filename (str): Name of file
        """

        fileutils.check_write_ok(self.opts.get('output_filename'),
                                 self.out_fmt,
                                 overwrite=self.opts.get('overwrite'))
        final_id = self._build()

        print("Writing system to %s with %d atoms comprising:\n"
              "  %d lipid molecules\n"
              "  %d water molecules\n"
              % (self.opts.get('output_filename'),
                 molutils.num_atoms_remaining(molid=final_id),
                 molutils.num_lipids_remaining(final_id, self.opts.get('lipid_sel')),
                 molutils.num_waters_remaining(molid=final_id)))
        fileutils.write_final_system(out_fmt=self.out_fmt,
                                     out_name=self.opts.get('output_filename'),
                                     molid=final_id, tmp_dir=self.tmp_dir,
                                     forcefield=self.opts.get('forcefield'),
                                     extra_topos=self.opts.get('extra_topos'),
                                     extra_params=self.opts.get('extra_params'),
                                     extra_streams=self.opts.get('extra_streams'),
                                     hmassrepartition=self.opts.get('hmassrepartition'),
                                     verbose=self.opts.get('debug_verbose')
                                    )
        molecule.delete(final_id)

    #==========================================================================

    def add_molecule(self, filename, desc):
        """
        Adds a molecule file to the system.

        Args:
          filename (str): File to load
          desc (str): Type of molecule: solute, solvent, etc

        Returns:
          True if the molecule was successfully added
        """

        molid = fileutils.load_solute(filename, tmp_dir=self.tmp_dir)
        self.molids[desc] = molid
        atomsel('all', molid=molid).beta = 1
        return True

    #==========================================================================

    def remove_molecule(self, desc):
        """
        Removes a molecule file from the system.

        Args:
          desc (str): Key for molecule to remove

        Returns:
          true if a molecule was deleted, false otherwise

        Raises:
          KeyError if there is no molecule with that key string
        """

        molid = self.molids[desc]
        try:
            molecule.delete(molid)
            return True
        except ValueError:
            return False
        return False

    #==========================================================================

    def convert_ions(self, salt_conc, molid):
        """
        Calculates the charge of the molecule and adds salt ions to get the
        desired concentration by converting water molecules to salt

        Args:
          salt_conc (float): Desired salt concentration in M
          cation (str): Cation to add
          anion (str): Anion to add
          molid (int): VMD molecule id to consider

        Returns:
          (int) number of ions added

        Raises:
          ValueError if invalid cation is specified
        """
        anion = self.opts["anion"]
        cation = self.opts["cation"]

        # Give existing cations correct nomenclature
        # TODO: I think this is where the bug with ions getting renamed is
        # molutils.set_cations(molid, cation)

        # Calculate number of salt ions needed
        # TODO: Handle non-1 charges here
        pos_ions_needed, neg_ions_needed, num_wat, total_cations, total_anions, \
        cation_conc, anion_conc = molutils.get_num_salt_ions_needed(molid,
                                                                    salt_conc,
                                                                    cation=cation,
                                                                    anion=anion)

        print("Solvent will be %d waters, %d %s (%.3f M), %d %s (%.3f M)" %
              (num_wat, total_cations, cation, cation_conc,
               total_anions, anion, anion_conc))

        print("Converting %d waters to %d %s ions and %d %s ions..." %
              (pos_ions_needed + neg_ions_needed,
               pos_ions_needed, cation, neg_ions_needed, anion))

        # Add the ions
        for _ in range(pos_ions_needed):
            add_salt_ion(cation, molid)
        for _ in range(neg_ions_needed):
            add_salt_ion(anion, molid)

        return pos_ions_needed + neg_ions_needed

    #==========================================================================

    def get_cell_size(self,
                      mem_buf, wat_buf,
                      molid=None,
                      filename=None,
                      zh_mem_full=_MEMBRANE_FULL_THICKNESS / 2.0,
                      zh_mem_hyd=_MEMBRANE_HYDROPHOBIC_THICKNESS / 2.0):
        """
        Gets the cell size of the final system given initial system and
        buffers. Detects whether or not a membrane is present. Sets the
        size of the system.

        Args:
          mem_buf (float) : Membrane (xy) buffer amount
          wat_buf (float) : Water (z) buffer amount
          molid (int) : VMD molecule ID to consider (can't use with filename)
          filename (str) : Filename of system to consider (can't use w molid)
          zh_mem_full (float) : Membrane thickness
          zh_mem_hyd (float) : Membrane hydrophobic region thickness

        Returns:
        return dx_sol, dy_sol, dx_tm, dy_tm, dz_full
          (float tuple): x solute dimension, y solute dimension,
            TM x solute dimension, TM y solute dimension, solute z dimension

        Raises:
          ValueError: if filename and molid are both specified
        """

        # Sanity check
        if filename is not None and molid is not None:
            raise ValueError("Specified molid and filename to get_cell_size")

        if filename is not None:
            top = molecule.get_top()
            molid = molecule.read(-1, 'mae', filename)
        elif molid is None:
            molid = molecule.get_top()

        # Some options different for water-only systems (no lipid)
        if self.water_only:
            solute_z = atomsel(self.solute_sel, molid=molid).z
            dx_tm = 0.0
            dy_tm = 0.0
            sol_solute = atomsel(self.solute_sel, molid)
        else:
            solute_z = atomsel(self.solute_sel, molid=molid).z
            tm_solute = atomsel('(%s) and z > %f and z < %f' % (self.solute_sel,
                                                                -zh_mem_hyd,
                                                                zh_mem_hyd), molid)
            if tm_solute:
                dx_tm = max(tm_solute.x)- min(tm_solute.x)
                dy_tm = max(tm_solute.y) - min(tm_solute.y)
            else:
                dx_tm = dy_tm = 0

            sol_solute = atomsel('(%s) and (z < %f or z > %f)' % (self.solute_sel,
                                                                  -zh_mem_hyd,
                                                                  zh_mem_hyd), molid)

        # Solvent invariant options
        dx_sol = max(sol_solute.x) - min(sol_solute.x)
        dy_sol = max(sol_solute.y) - min(sol_solute.y)

        if self.opts.get('user_x'):
            self.size[0] = self.opts['user_x']
        else:
            self.size[0] = max(dx_tm + 2.*mem_buf, dx_sol + 2.*wat_buf)
        if self.opts.get('user_y'):
            self.size[1] = self.opts['user_y']
        else:
            self.size[1] = max(dy_tm + 2.*mem_buf, dy_sol + 2.*wat_buf)

        # Z dimension. If there's a membrane, need to account for asymmetry
        # in the Z dimension where the protein could be uneven in the membrane
        # or even peripheral
        if self.opts.get('user_z'):
            self.size[2] = self.opts['user_z']
            buf = (self.opts['user_z'] - max(solute_z) + min(solute_z))/2
            self._zmax = max(solute_z) + buf
            self._zmin = min(solute_z) - buf
            if zh_mem_full > self._zmax or -zh_mem_full < self._zmin:
                raise DabbleError("Specified user z of %f is too small to "
                                  "accomodate protein and membrane!"
                                  % self.opts['user_z'])
        else:
            if self.water_only:
                self._zmax = max(solute_z) + wat_buf
                self._zmin = min(solute_z) - wat_buf
            else:
                self._zmax = max(max(solute_z)+wat_buf, zh_mem_full)
                self._zmin = min(min(solute_z)-wat_buf, -zh_mem_full)
            self.size[2] = self._zmax - self._zmin

        # Cleanup temporary file, if read in
        if filename is not None:
            molecule.delete(molid)
            if top != -1:
                molecule.set_top(top)

        return dx_sol, dy_sol, dx_tm, dy_tm, max(solute_z)-min(solute_z)

    #==========================================================================
    #                            Private methods                              #
    #==========================================================================

    def _set_solute_sel(self, molid):
        """
        Gets the list of resids uniquely corresponding to the in the system
        and sets the value of the solute_sel attribute.
        Does this separately by chain since sometimes resids can be the
        same across different chains. This selection can be used to pull out
        the solute once other things are added later.
        This assumes that the solute is the only thing in the system right now.

        Args:
          molid (int) : VMD molecule ID to get the selection from

        Returns:
          (str) : VMD atom selection for these residues
        """
        # Temporary fix for chain W in input file
        if atomsel('chain W'):
            print("WARNING: Renaming crystal water chain to X, temporary bugfix")
            atomsel('chain W').chain = 'X'

        chains = set(atomsel('all', molid=molid).chain)
        sel = ""
        while chains:
            chn = chains.pop() # have to handle first separately because of or
            sel += "(chain %s and resid " % (chn) + \
                   " ".join(["'%d'" % i for i in \
                       set(atomsel('chain %s' % chn, molid=molid).resid)]) + \
                   ")"
            if chains:
                sel += " or "

        self.solute_sel = sel
        return sel

        #==========================================================================

    def _set_cell_to_square_prism(self, molid):
        """
        Sets the periodic box to be a square prism of specified dimension

        Args:
          molid (int): VMD molecule ID to consider
        """
        old_top = molecule.get_top()
        molecule.set_top(molid)
        molecule.set_periodic(-1, -1,
                              self.size[0], self.size[1], self.size[2],
                              90.0, 90.0, 90.0)
        molecule.set_top(old_top)

    #==========================================================================

    def _remove_z_residues(self, molid):
        """
        Removes residues in the +-Z direction in the system. Used to chop off
        extra waters away from the protein to keep the system the desired
        size. Changes apparent on next write.

        Args:
          molid (int): VMD molecule id to use

        Returns:
          (int) number of atoms deleted
        """
        return _remove_residues('(not (%s)) and noh and abs(z) > %f'
                                % (self.solute_sel, self.size[2] / 2.0),
                                molid=molid)

    #==========================================================================

    def _add_water(self, molid):
        """
        Adds water residues in the +- Z direction in the system. Used if there
        are not enough waters in the initial membrane system. Creates a new
        molecule with waters added.

        Args:
            molid (int): VMD molecule id to use

        Returns:
            (int) VMD molecule id with water added

        Raises:
            ValueError if water buffer is None
        """

        if not self.opts.get('wat_buffer'):
            raise ValueError("Water buffer undefined")

        # Check the +Z direction
        zup = self.opts['wat_buffer'] + \
              max(atomsel(self.solute_sel).z) - \
              max(atomsel("not (%s or %s)" % (self.solute_sel, self.opts['lipid_sel'])).z)
        # Check the -Z direction
        zdo = self.opts['wat_buffer'] + \
              min(atomsel("not (%s or %s)" % (self.solute_sel, self.opts['lipid_sel'])).z) \
              - min(atomsel(self.solute_sel).z)

        # Load water
        wat_path = self.opts['membrane_system'] = resource_filename(__name__, \
                "lipid_membranes/tip3pbox.mae")
        self.add_molecule(wat_path, 'water')
        to_combine = [molid]

        # Handle adding water above the protein
        # When moving add a 0.5 less so there isn't a gap
        if zup > 0:
            print("Adding %f A water above the solute..." % zup)
            self.molids['wtmp'], _ = \
                tile_membrane_patch(self.molids['water'],
                                    [self.size[0], self.size[1], zup],
                                    self.tmp_dir, allow_z_tile=True)

            move = max(atomsel("not (%s)" % self.solute_sel, molid=molid).z) - \
                    min(atomsel(molid=self.molids['wtmp']).z) - 0.5
            atomsel(molid=self.molids['wtmp']).moveby((0, 0, move))
            self.molids['wats_up'] = molutils.center_system(molid=self.molids['wtmp'],
                                                            tmp_dir=self.tmp_dir,
                                                            center_z=False)
            self.remove_molecule('wtmp')
            to_combine.append(self.molids['wats_up'])

        # Handle adding water below the protein
        # When moving add a 0.5 less so there isn't a gap
        if zdo > 0:
            print("Adding %f A water below the solute..." % zdo)
            self.molids['wtmp'], _ = \
                tile_membrane_patch(self.molids['water'],
                                    [self.size[0], self.size[1], zdo],
                                    self.tmp_dir, allow_z_tile=True)

            move = min(atomsel("not (%s)" % self.solute_sel, molid=molid).z) - \
                    max(atomsel(molid=self.molids['wtmp']).z) + 0.5
            atomsel(molid=self.molids['wtmp']).moveby((0, 0, move))
            self.molids['wats_down'] = molutils.center_system(molid=self.molids['wtmp'],
                                                              tmp_dir=self.tmp_dir,
                                                              center_z=False)
            self.remove_molecule('wtmp')
            to_combine.append(self.molids['wats_down'])

        # Combine and return a new molecule
        newid = molutils.combine_molecules(input_ids=to_combine,
                                           tmp_dir=self.tmp_dir)
        self.remove_molecule('water')
        return newid

    #==========================================================================

    def _trim_water(self, molid):
        """
        Removes water residues in the +- Z direction in the system. Used
        to chop off extra waters from the protein that are past a certain cutoff
        distance. Changes are apparent on next write.

        Args:
          molid (int): VMD molecule id to use

        Returns:
          int Number of atoms deleted

        Raises:
          ValueError if water buffer is None
        """

        if not self.opts.get('wat_buffer'):
            raise ValueError("Water buffer undefined")

        # Remove waters in the Z direction
        # Check if we are trimming to absolute size and set z buf if so
        total = 0
        total = _remove_residues('(not (%s) and not (%s)) and noh and z > %f' % \
                                 (self.solute_sel, self.opts['lipid_sel'],
                                  self._zmax), molid=molid)
        total += _remove_residues('(not (%s) and not (%s)) and noh and z < %f' %
                                  (self.solute_sel, self.opts['lipid_sel'],
                                   self._zmin), molid=molid)

        # Trim in the XY direction if it's a pure water system
        # lipid trimming is done in trim_xy_residues and takes into account
        # lipid center, etc.
        if self.water_only:
            xcoord = atomsel(self.solute_sel).x
            buf = (self.size[0] - max(xcoord) + min(xcoord))/2.
            total += _remove_residues('(not (%s)) and noh and x > %f' %
                                      (self.solute_sel, max(xcoord) + buf),
                                      molid=molid)
            total += _remove_residues('(not (%s)) and noh and x < %f' %
                                      (self.solute_sel, min(xcoord) - buf),
                                      molid=molid)

            ycoord = atomsel(self.solute_sel).y
            buf = (self.size[1] - max(ycoord) + min(ycoord))/2.
            total += _remove_residues('(not (%s)) and noh and y > %f' %
                                      (self.solute_sel, max(ycoord) + buf),
                                      molid=molid)
            total += _remove_residues('(not (%s)) and noh and y < %f' %
                                      (self.solute_sel, min(ycoord) - buf),
                                      molid=molid)

        return total

    #==========================================================================

    def _remove_xy_lipids(self, molid):
        """
        Removes residues in the +-XY direction in the system. Used to chop off
        lipids that are protruding outside of the box dimensions.

        Args:
          molid (int): VMD molecule id to use

        Returns:
          (int) number of atoms deleted

        Raises:
          ValueError if the 'lipid' only contains hydrogens
        """

        # Select residues that are outside the box
        half_x_size = self.size[0] / 2.0
        half_y_size = self.size[1] / 2.0
        box_sel_str = 'abs(x) > %f or abs(y) > %f' % (half_x_size, half_y_size)

        # Identify lipids that have some part outside of the box
        suspicious_lipid_residues = list(set(atomsel('(%s) and (%s)' % \
               (self.opts['lipid_sel'], box_sel_str), molid=molid).residue))
        bad_lipids = []

        # Delete lipids whose center is too far out of the box, keep others
        for i in suspicious_lipid_residues:
            lipid_center = atomsel('noh and residue %s' % str(i),
                                   molid=molid).center()
            # Sanity check
            if not lipid_center:
                raise DabbleError("No heavy atoms found in suspicious residue %s"
                                  "\nCheck your input file." % str(i))

            if abs(lipid_center[0]) > half_x_size or \
               abs(lipid_center[1]) > half_y_size:
                bad_lipids.append(i)
        lipid_headgroup_sel = 'residue ' + ' '.join([str(l) for l in bad_lipids])

        # Do the deletion
        removal_sel_str = '(%s) or not (%s)' % (lipid_headgroup_sel,
                                                self.opts['lipid_sel'])
        total = _remove_residues('noh and (%s) and (%s) and not (%s)' %
                                 (box_sel_str,
                                  removal_sel_str,
                                  self.solute_sel),
                                 molid=molid)
        return total

    #==========================================================================

    def _remove_overlapping_residues(self,
                                     lipid_sel,
                                     molid,
                                     lipid_friendly_sel=None,
                                     dist=1.75):
        """
        Removes residues that are overlapping. For example, when the protein
        system is combined with the membrane, remove the lipids that are now in
        the same place as the protein. Changes will be written on next call
        to write to file.

        Args:
          lipid_sel (str): VMD atom selection for the lipids
          molid (int): VMD molecule to remove from
          lipid_friendly_sel (str): VMD atom selection for lipids that are
            allowed to be much closer to the protein, or None
          dist (float): Minimum distance between atoms, defaults to 1.75 A

        Returns:
          (int) number of atoms removed due to clashes
        """
        # Select and remove solvent molecules that are clashing
        # If water only, remove waters whose oxygen is within 2.5A of a heavy
        # atom in the solute (doesn't grab crystal waters in original structure!)
        if self.water_only:
            clashing_sel = "noh and (not (%s)) and " \
                           "(pbwithin %f of (noh and (%s)))" \
                           % (self.solute_sel, max(dist, 2.5), self.solute_sel)
        else:
            clashing_sel = 'not (%s) and noh and not (%s) and ' \
            '(pbwithin %f of (noh and (%s)))' % (lipid_sel, self.solute_sel,
                                                 dist, self.solute_sel)
        total = _remove_residues(clashing_sel, molid=molid)

        # Select and remove lipid molecules that are clashing
        if lipid_friendly_sel is not None:
            clashing_sel_lipid = '(%s) and noh and not (%s) and ' \
            'pbwithin %f of (noh and (%s) and not (%s))' \
                    % (lipid_sel, self.solute_sel, dist, self.solute_sel,
                       lipid_friendly_sel)
        else:
            clashing_sel_lipid = '(%s) and noh and not (%s) and ' \
            'pbwithin %f of (noh and (%s))' % (lipid_sel, self.solute_sel,
                                               dist, self.solute_sel)
        total += _remove_residues(clashing_sel_lipid, molid=molid)
        return total

    #==========================================================================

    def _remove_lipids_near_rings(self,
                                  lipid_sel,
                                  molid,
                                  ring_sel='noh and resname HID HIE HIP '
                                           'HIS PHE TRP TYR and not backbone',
                                  dist=1.75):
        """
        Deletes lipids that could stick through protein rings. Changes
        will be written on the next call to write to a file

        Args:
          lipid_sel (str): VMD atom selection for the lipids
          molid (int): VMD molecule id to look at
          ring_sel (str): VMD atom selection for ring amino acids to check
          dist (float): Minimum distance between atoms, defaults to 1.75 A

        Returns:
          (int) number of atoms removed due to sticking through rings
        """

        solute_ring_sel = '(%s) and (%s)' % (self.solute_sel, ring_sel)
        return _remove_residues('noh and (%s) and not (%s) and '
                                'pbwithin %f of (noh and (%s))'
                                % (lipid_sel, self.solute_sel,
                                   dist, solute_ring_sel),
                                molid=molid)

    #==========================================================================

    def _remove_lipid_boundary_clash(self, pointy_type, ring_type, # pylint: disable=no-self-use
                                     molid, dist=1.0):
        """
        Deletes lipids that are too close to other lipids at the periodic
        boundary

        Args:
          pointy_type (str): VMD atom selection for pointy type
          ring_type (str): VMD atom selection for ring type
          molid (int): VMD molecule id to look at
          dist (float): Minimum distance between atoms, defaults to 1.0 A

        Returns:
          (int) number of atoms removed due to boundary clash
        """
        # TODO: arent these selections just the same?
        sel1 = 'noh and (%s) and pbwithin %f of noh and (%s)' % (pointy_type,
                                                                 dist,
                                                                 ring_type)
        sel2 = 'noh and (%s) and pbwithin %f of noh and (%s)' % (ring_type,
                                                                 dist,
                                                                 pointy_type)
        return _remove_residues('(%s) or (%s)' % (sel1, sel2), molid=molid)

    #==========================================================================

    def _remove_clashing_lipids(self, molid, lipid_sel, lipid_friendly_sel):
        """
        Removes all types of clashing lipids in a given molecule

        Args:
          molid (int): VMD molecule id to look at
          lipid_sel (str): VMD atom selection for lipids

        Returns:
          (int) total number of atoms removed
        """

        # Remove lipids stuck through aromatic rings
        ring_lips = self._remove_lipids_near_rings(lipid_sel,
                                                   molid=molid)
        # TODO
        clash_lips = 0
        protein_lips = 0
        solute_lips = 0
        if self.opts.get('clash_lipids'):
            x_edge_dim = self.size[0]*0.5*0.9
            y_edge_dim = self.size[1]*0.5*0.9
            pointy_lipid_sel = '((%s) and not (%s)) and (abs(x) > %f or ' \
            'abs(y) > %f)' % (lipid_sel, self.opts['clash_lipids'],
                              x_edge_dim, y_edge_dim)
            ring_lipid_sel = '%s and (abs(x) > %f or abs(y) > %f)' \
            % (self.opts['clash_lipids'], x_edge_dim, y_edge_dim)
            clash_lips = self._remove_lipid_boundary_clash(pointy_lipid_sel,
                                                           ring_lipid_sel,
                                                           molid)

            #TODO replace with its own thing not near rings
            protein_lips = self._remove_lipids_near_rings(self.opts['clash_lipids'], \
                ring_sel='noh and protein and not backbone', dist=1.25, \
                molid=molid)

        solute_lips = self._remove_xy_lipids(molid)

        print("Removed %d atoms" % (solute_lips + ring_lips + clash_lips +
                                    protein_lips))
        print("  %4d atoms from lipids running into the solute\n"
              "  %4d atoms from lipids going through HIS, PHE, TYR, TRP\n"
              "  %4d atoms from one lipid going through another\n"
              "  %4d atoms from lipids getting stuck in a sidechain"
              % (solute_lips, ring_lips, clash_lips, protein_lips))

#==========================================================================

    def _orient_solute(self, molid):
        """
        Orients the solute. Can either move it explicitly in the z direction,
        or align to an OPM structure.

        Args:
          molid (int): VMD molecule ID to orient
          z_move (float): Amount to move in the Z direction
          z_rotation (float): Amount to rotate membrane relative to protein,
            can just take this straight from the OPM website value
          opm_pdb (str): Filename of OPM structure to align to
          opm_align (str): Atom selection string to align
          tmp_dir (str): Directory to put temporary files in

        Returns:
          (int) VMD molecule ID of oriented system

        Raises:
          ValueError if movement and alignment arguments are both specified
        """


        # Check that OPM and alignment aren't both specified
        if self.opts.get('opm_pdb') and \
                (self.opts.get('z_move') != 0 or self.opts.get('z_rotation') != 0):
            raise DabbleError("ERROR: Cannot specify an OPM pdb and "
                              "manual orientation information")

        if self.opts.get('opm_pdb'):
            opm = molecule.load('pdb', self.opts['opm_pdb'])
            moveby = atomsel('protein and backbone', molid=molid).fit( \
                             atomsel(self.opts.get('opm_align'), molid=opm))
            atomsel('all', molid=molid).move(moveby)
            molecule.delete(opm)
            return molid

        if self.opts.get('z_move'):
            atomsel('all', molid=molid).moveby((0, 0, self.opts['z_move']))
            if not self.opts.get('z_rotation'):
                return molid

        if self.opts.get('z_rotation'):
            trans.resetview(molid) # View affect rotation matrix, now it's I
            # This is negative because we want membrane flat along the z-axis,
            # and OPM lists the membrane rotation relative to the protein
            theta = math.radians(-1*self.opts['z_rotation'])
            # Rotation matrix in row order with 4th dimension just from I
            # pylint: disable=bad-whitespace, bad-continuation
            rotmat = [ math.cos(theta), -1*math.sin(theta), 0, 0,
                       math.sin(theta),   math.cos(theta),  0, 0,
                             0        ,         0,          1, 0,
                             0        ,         0,          0, 1 ]
            # pylint: enable=bad-whitespace, bad-continuation
            trans.set_rotation(molid, rotmat)
            return molid

        # Center the system according to VMD's internal metric, then
        # move the protein in the xy plane so that there is equal padding
        # on either side
        molid = molutils.center_system(molid=molid, tmp_dir=self.opts.get('tmp_dir'),
                                       center_z=self.water_only)
        system = atomsel('all', molid=molid)
        tx = (-max(system.x) - min(system.x))/2.
        ty = (-max(system.y) - min(system.y))/2.
        temp_mae = tempfile.mkstemp(suffix='.mae',
                                    prefix='dabble_centered',
                                    dir=self.opts.get('tmp_dir'))[1]
        system.moveby((tx, ty, 0))
        system.write('mae', temp_mae)
        molecule.delete(molid)
        new_id = molecule.load('mae', temp_mae)
        return new_id

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                           MODULE FUNCTIONS
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _find_convertible_water_molecule(molid, # pylint: disable=invalid-name
                                     water_sel='resname HOH WAT TIP3',
                                     min_ion_dist=5.0):
    """
    Finds a water molecule that can be converted to an ion

    Args:
      molid (int): VMD molid to look at
      water_sel (str): VMD atom selection for water
      min_ion_dist (float): Minimum distance between ionds

    Returns:
      (int) Atom index of a water oxygen that is convertible

    Raises:
      ValueError if no convertible water molecules are found
    """

    inclusion_sel = 'beta 1 and noh and (%s)' % water_sel
    exclusion_sel = 'beta 1 and not (%s)' % water_sel
    sel = atomsel('(%s) and not pbwithin %f of (%s)' \
                  % (inclusion_sel, min_ion_dist, exclusion_sel), molid).index
    if not sel:
        raise DabbleError("No convertible water molecules found in %s" % sel)

    return sel[random.randint(0, len(sel)-1)]

#==========================================================================

def _convert_water_molecule_to_ion(molid, atom_id, element):
    """
    Converts a water molecule to an ion, deleting the hydgrogens and
    placing the ion where the water oxygen was

    Args:
      molid (int): VMD molecule to operate on
      atom_id (int): Atom index of water oxygen to change to ion
      element (str): Ion to apply

    Raises:
      ValueError: if invalid element specified
      ValueError: if atom id is not of an oxygen water
    """

    element_received = atomsel('index %d' % atom_id).element[0]
    if element_received != 'O':
        raise ValueError("Received non-water oxygen to convert, id %d "
                         "element %s" % (atom_id, element_received))

    molutils.set_ion(molid, atom_id, element)
    _remove_atoms('element H and same residue as index %d' % atom_id, molid)

#==========================================================================

def _remove_atoms(sel, molid):
    """
    Marks specified atoms for removal. IMPORTANT - atoms are not actually
    deleted. Changes will be written on next call to write!

    Args:
      sel (str): VMD atom selection string to remove
      molid (int): VMD molecule ID to consider

    Returns:
      (int): The number of atoms removed
    """

    marked = atomsel('beta 1 and (%s)' % sel, molid=molid)
    marked.beta = 0
    return len(marked)

#==========================================================================

def _remove_residues(sel, molid):
    """
    Marks all residues containing an atom in the selection for deletion.
    IMPORTANT - atoms are not actually deleted until next call to write!

    Args:
      sel (str): VMD atom selection string to remove
      molid (int): VMD molecule ID to consider

    Returns:
      (int): The number of atoms removed
    """

    return _remove_atoms('same residue as (%s)' % sel, molid)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                            PUBLIC FUNCTIONS                             #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def add_salt_ion(element, molid):
    """
    Changes a water molecule to a salt ion.

    Args:
      element (str): ion to add
      molid (int): VMD molecule id to consider

    Raises:
      AssertionException : if element is not supported

    Returns:
      (int) the index of the water molecule that was replaced
    """
    atom_id = _find_convertible_water_molecule(molid)
    _convert_water_molecule_to_ion(molid, atom_id, element)
    return atom_id

#==========================================================================

def tile_membrane_patch(input_id, min_size, tmp_dir, allow_z_tile):
    """
    Tiles a system in the x and y dimension (z currently unsupported) to
    make a larger system.

    Args:
      input_id (int): VMD molecule id to tile
      min_size (array of 3 floats): Final system X, Y, Z dimension
      tmp_dir (str): Directory in which to put tiled molecule
      allow_z_tile (bool): Whether to allow tiling in the Z direction

    Returns:
      (int) output molecule id
      (int 3x) number of times tiled in x, y, z direction
    """

    mem_dimensions = np.array(molutils.get_system_dimensions(molid=input_id))
    times_x, times_y, times_z = [int(times) for times in \
            np.ceil(min_size / mem_dimensions)]

    # Disallow tiling in Z direction
    if not allow_z_tile:
        times_z = 1

    # If there is not enough water in the Z direction, it will be added later

    if times_x == 1 and times_y == 1 and times_z == 1:
        output_id = input_id
    else:
        output_id = molutils.tile_system(input_id, times_x, times_y, times_z,
                                         tmp_dir=tmp_dir)

    return output_id, (times_x, times_y, times_z)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
