"""
This module contains the CharmmWriter class and associated methods,
which outputs a psf/pdb file with CHARMM names and parameters.
It does this by converting atom names to CHARMM names, writing
intermediate files as necessary to invoke the vmd psfgen plugin.

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
import sys
import os
import tempfile
from pkg_resources import resource_filename

import vmd
import molecule
from atomsel import atomsel

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                CONSTANTS                                    #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_acids = ('ACE ALA ARG ASN ASP CYS CYX GLN GLU GLY HIE HIS HSP HSE '
          'HSD ILE LEU LYS MET NMA PHE PRO SER THR TRP TYR VAL')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                   CLASSES                                   #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CharmmWriter(object):
    """
    An object that handles all the conversions to a psf file
    by interfacing with psfgen.

    Writes a pdb/psf file pair from the current molecule using the
    CHARMM36 topology and atom names/types. Interfaces with psfgen by
    dynamically generating the .tcl file that psfgen takes as input.
    Prompts the user for additional topology files and helps with
    matching atom names that cannot be automatically translated to the
    charmm naming conventions.

    Attributes:
      file (file handle): Temporary file to write TCL script that invokes
          psfgen
      tmp_dir (string): Directory where temporary files are stored
      psf_name (str): Prefix for the pdb/psf output files, extension will
          be appended
      molid (str,optional): the VMD molecule id to write. Defaults to 0.
      lipid_sel (str,optional): atomselect string describing what should count
          as "lipid".  Defaults to "lipid"
      topologies (list of str): Topology files that were used in creating the
          psf
      prompt_topos (bool): Whether to ask for more topology files

    """

    #==========================================================================

    def __init__(self, tmp_dir, molid, lipid_sel="lipid", extra_topos=None):
        # Create TCL temp file and directory
        self.tmp_dir = tmp_dir 
        self.filename = tempfile.mkstemp(suffix='.tcl', prefix='dabble_psfgen',
                                         dir=self.tmp_dir)[1]
        self.lipid_sel = lipid_sel
        self.file = open(self.filename, 'w')
        self.molid = molid
        self.psf_name = ""
        # Default parameter sets
        self.topologies = [
            resource_filename(__name__, "charmm_parameters/top_all36_caps.rtf"),
            resource_filename(__name__, "charmm_parameters/top_water_ions.rtf"),
            resource_filename(__name__, "charmm_parameters/top_all36_cgenff.rtf"),
            resource_filename(__name__, "charmm_parameters/top_all36_prot.rtf"),
            resource_filename(__name__, "charmm_parameters/top_all36_lipid.rtf"),
            resource_filename(__name__, "charmm_parameters/top_all36_carb.rtf")
            ]
        if extra_topos is not None:
            self.topologies.extend(extra_topos)
            self.prompt_topos = False
        else:
            self.prompt_topos = True

    #=========================================================================

    def write(self, psf_name):
        """
        Writes the pdb/psf file.

        Args:
          psf_name (str): Prefix for the pdb/psf output files, extension
            will be appended

        Returns:
          topologies (list of str): Topology files that were used in creating
              the psf
        """
        # Clean up all temp files from previous runs if present
        # An earlier check will exit if it's not okay to overwrite here
        self.psf_name = psf_name
        try:
            os.remove('%s.pdb'% self.psf_name)
            os.remove('%s.psf'% self.psf_name)
        except OSError:
            pass

        #  Finds the psfgen package and sets the output file name
        string = '''
        set dir [file join $env(VMDDIR) plugins [vmdinfo arch] tcl psfgen1.6]
        package ifneeded psfgen 1.6.2 [list load [file join $dir libpsfgen.so]]
        package require psfgen
        set output "%s"
        resetpsf
        ''' % self.psf_name
        self.file.write(string)

        # Put our molecule on top
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        # Ask the user for additional topology files
        if self.prompt_topos:
            sys.stdout.flush()
            print("\nCurrently using the following topology files:")
            for top in self.topologies:
                print("  - %s" % top.split("/")[-1])

            print("Enter the path to the filename(s) from the current working "
                  "directory, separated by a comma, of any additional rtf files "
                  "you wish to use.\n")
            sys.stdout.flush()
            inp = raw_input('> ')
            if inp:
                self.topologies.extend(inp.split(','))
        else:
            print("Using the following topologies:")
            for top in self.topologies:
                print("  - %s" % top.split("/")[-1])

        self.file.write('\n')
        for top in self.topologies:
            self.file.write('   topology %s\n' % top)

        # Mark all atoms as unsaved with the user field
        atomsel('all', molid=self.molid).set('user', 1.0)

        # Now ions if present, changing the atom names
        if len(atomsel('element Na Cl K', molid=self.molid)) > 0:
            self._write_ion_blocks()

        # Save water 10k molecules at a time
        if len(atomsel('water', molid=self.molid)):
            self._write_water_blocks()

        # Now lipid
        if len(atomsel(self.lipid_sel)):
            self._write_lipid_blocks()

        if not len(atomsel('resname %s' % _acids, molid=self.molid)):
            print("\nINFO: Didn't find any protein.\n")

        # Pull out the protein, one fragment at a time
        for frag in set(atomsel('resname %s' % _acids).get('fragment')):
            temp = tempfile.mkstemp(suffix='_P%s.mae' % frag,
                                    prefix='psf_prot_', dir=self.tmp_dir)[1]
            fragment = atomsel('fragment %s' % frag, molid=self.molid)
            fragment.write('mae', temp)
            prot_molid = molecule.load('mae', temp)

            self._write_protein_blocks(seg='P%s'%frag,
                                       prot_molid=prot_molid)
            molecule.delete(prot_molid)
            fragment.set('user', 0.0)

        # Detect disulfide bridges and add appropriate patch lines
        self._set_disulfide_bridges()
        # End protein

        # Check if there is anything else and let the user know about it
        leftovers = atomsel('user 1.0', molid=self.molid)
        if len(leftovers) > 0:
            print("Found extra ligands: %s" % set(leftovers.get('resname')))
        for lid in set(leftovers.get('residue')):
            self._write_ligand_blocks(lid)

        # Write the output files and run
        string = '''
        writepsf x-plor cmap ${output}.psf
        writepdb ${output}.pdb'''
        self.file.write(string)
        self.file.close()

        from VMD import evaltcl
        evaltcl('play %s' % self.filename)
        self._check_psf_output()

        # Reset top molecule
        molecule.set_top(old_top)

        return self.topologies

    #=========================================================================
    #                           Private methods                              #
    #=========================================================================

    def _write_protein_blocks(self, seg, prot_molid):
        """
        Writes a temporary protein PDB file with correct atom names for psfgen,
        in atom order.

        Args:
          seg (str): VMD segment to write
          protmolid (int): VMD molecule ID to write, usually a temp protein
        """

        # Put molid on top to simplify atom selections
        old_top = molecule.get_top()
        molecule.set_top(prot_molid)
        print("Writing protein file\n")

    # TODO There is a bug where somewhere the NMA gets named ASP
        #atomsel().write('mae','duh.mae')
        #quit()

        # Renumber residues starting from 1
        atomsel('all').set('user', 1.0)
#       residues = set(atomsel('all').get('residue'))
#       resnum = 1
#       while len(residues):
#           atomsel('residue %s'% residues.pop()).set('resid', resnum)
#           resnum += 1

        # Check the protein has hydrogens
        if len(atomsel('element H')) == 0:
            print("\n\nERROR: There are no hydrogens on your protein")
            print("       You need to add them before parameterizing "
                  "because psfgen cannot be trusted to do it correctly\n")
            quit(1)

        # Terminal residue ACE
        ace_names = {'CH3' :'CAY', 'O'   :'OY',
                     'HH31':'HY1', 'HH32':'HY2', 'HH33':'HY3',
                     'H1'  :'HY1', 'H2'  :'HY2', 'H3'  :'HY3',
                     '1H'  :'HY1', '2H'  :'HY2', '3H'  :'HY3'}
        for name in ace_names:
            atomsel('resname ACE and name %s'%name).set('name', ace_names[name])

        # Terminal residue NMA
        nma_names = {'HN'  :'HNT', 'H'   :'HNT',
                     'CT3' :'CAT', 'CA'  :'CAT', 'CH3' :'CAT',
                     'HA1' :'HT1', 'HA2' :'HT2', 'HA3' :'HT3',
                     '1HA' :'HT1', '2HA' :'HT2', '3HA' :'HT3',
                     'HH31':'HT1', 'HH32':'HT2', 'HH33':'HT3'}
        for name in nma_names:
            atomsel('resname NMA and name %s'%name).set('name', nma_names[name])


        # Determine protonation states of those called HIS based on atom names
        atomsel('resname HIS and same residue as name HE2').set('resname',
                                                                'HSE')
        atomsel('resname HIS and same residue as name HD1').set('resname',
                                                                'HSD')
        # NOTE: This atomsel MUST come after the previous two!
        atomsel('resname HIS and same residue as name HE2 '
                'and same residue as name HD1').set('resname', 'HSP')

        # Isoleucine
        iso_names = {'CD1' :  'CD',
                     'HD11': 'HD1', 'HD12':'HD2', 'HD13':'HD3',
                     'HG12':'HG11', 'HG13':'HG12'}
        for name in iso_names:
            atomsel('resname ILE and name %s'%name).set('name', iso_names[name])

        # Lysine last atom
        atomsel('resname LYS and name HZ').set('name', 'HZ3')

        patches = ''

        # Histidine naming convention
        atomsel('resname HID').set('resname', 'HSD')
        atomsel('resname HIE').set('resname', 'HSE')
        atomsel('resname HIP').set('resname', 'HSP')

        # Glutamine check all residues for protonation
        # Can't use dictionary for names since two of them must be swapped
        # Must loop by resid, not residue, to get correct PATCH statement index
        glups = set(atomsel('resname GLU GLH GLUP').get('resid'))
        for resid in glups:
            atomsel('resid %s' % resid).set('resname', 'GLU')
            if "HE1" in atomsel('resid %s' % resid).get('name'):
                atomsel('resid %s and name HE1'% resid).set('name', 'HE2')
                atomsel('resid %s and name OE1'% resid).set('name', 'temp')
                atomsel('resid %s and name OE2'% resid).set('name', 'OE1')
                atomsel('resid %s and name temp'% resid).set('name', 'OE2')
                patches += 'patch GLUP %s:%d\n' % (seg, resid)
            elif "HE2" in atomsel('resid %s' % resid).get('name'):
                patches += 'patch GLUP %s:%d\n' % (seg, resid)
            elif "HXT" in atomsel('resid %s' % resid).get('name'):
                atomsel('resid %s and name HXT' % resid).set('name', 'HE2')
                patches += 'patch GLUP %s:%d\n' % (seg, resid)

        # Aspartate check each residue to see if it's protonated
        # Can't use dictionary for names since two of them must be swapped
        # Must loop by resid, not residue, to get correct PATCH statement index
        asps = set(atomsel('resname ASP ASH ASPP').get('resid'))
        for resid in asps:
            atomsel('resid %s'% resid).set('resname', 'ASP')
            if "HD1" in atomsel('resid %s' % resid).get('name'):
                atomsel('resid %s and name HD1' % resid).set('name', 'HD2')
                atomsel('resid %s and name OD1' % resid).set('name', 'temp')
                atomsel('resid %s and name OD2' % resid).set('name', 'OD1')
                atomsel('resid %s and name temp' % resid).set('name', 'OD2')
                patches += 'patch ASPP %s:%d\n' % (seg, resid)
            if "HD2" in atomsel('resid %s' % resid).get('name'):
                patches += 'patch ASPP %s:%d\n' % (seg, resid)
# ROBIN: disabled for now-- too many false positives
#        # Check if an N terminal patch is needed
#        if "ACE" not in atomsel().get('resname'):
#            resid = min(atomsel().get('resid'))
#            index = atomsel('name N and resid %s' % resid).get('index')[0]
#            v = self._get_bonded_atoms(prot_molid, index)
#            v.sort()
#            if cmp(v, ["C","H","H","H"]) == 0:
#                print("INFO: Found N-terminal resid %d" % resid)
#                patches += 'patch NTER %s:%d\n' % (seg, resid)
#
#        # Check if a C terminal patch is needed
#        if "NMA" not in atomsel().get('resname'):
#            resid = max(atomsel().get('resid'))
#            index = atomsel('name N and resid %s' % resid).get('index')[0]
#            v = self._get_bonded_atoms(prot_molid, index)
#            v.sort()
#            if cmp(v, ["C","O","N"]):
#                print("INFO: Found amidated C-terminal resid %d" % resid)
#                patches += 'patch CT2 %s:%d\n' % (seg, resid)
#            elif cmp(v, ["C","O","O"]):
#                print("INFO: Found C-terminal resid %d" % resid)
#                patches += 'patch CTER %s:%d\n' % (seg,resid)

        # Methionine hydrogen names
        atomsel('name H2 H1 and resname MET').set('name', 'HN')

        # Serine and cysteine hydrogen names
        atomsel('name HG and resname SER CYS').set('name', 'HG1')

        # Hydrogens can be somewhat residue-dependent, but only check Hs
        # on known amino acid residues so that the names on nonstandard amino
        # acids are never changed
        h_names = {'H'  : 'HN', 'H2' :'HN', 'H1' : 'HN',
                   'HA2':'HA1',
                   'HA3':'HA2', 'HG2':'HG1',
                   'HG3':'HG2'}
        for hyd in h_names:
            atomsel('resname %s and name %s'
                    % (_acids, hyd)).set('name', h_names[hyd])

        # These two statements must execute in this specific order
        atomsel('resname %s and name HB2 and '
                'not resname ALA'% _acids).set('name', 'HB1')
        atomsel('resname %s and name HB3 and '
                'not resname ALA'% _acids).set('name', 'HB2')

        # Some residue-specific stuff
        atomsel('resname %s and name HD2 and not (resname '
                'TYR ILE HSP HSE HSD ASP PHE)' % _acids).set('name', 'HD1')
        atomsel('resname %s and name HD3 and not (resname '
                'TYR ILE HIS HSE HSD ASP PHE)'% _acids).set('name', 'HD2')
        atomsel('resname %s and name HE2 and not (resname '
                'TYR MET TRP HSP HSE HSD GLU PHE)'% _acids).set('name', 'HE1')
        atomsel('resname %s and name HE3 and not (resname '
                'TYR MET TRP HSP HSE HSD PHE)'% _acids).set('name', 'HE2')

        # Final chance to fix unrecognized atom names for non-protein residues
        others = set(atomsel('not resname %s' % _acids).get('residue'))
        if len(others):
            print("WARNING: Found non-protein residues in protein...")

        for residue in others:
            while not self._find_residue_in_rtf(residue=residue, molid=prot_molid):
                res = atomsel('residue %s' % residue).get('resname')[0]
                print("\nERROR: Residue name %s wasn't found in any input "
                      "topology. Would you like to rename it?\n" % res)
                sys.stdout.flush()
                newname = raw_input("New residue name or CTRL+D to quit > ")
                sys.stdout.flush()
                atomsel('residue %s' % residue).set('resname', newname)

        # If ACE and NMA aren't present, prompt for the residue name of the patch to 
        # apply, since auto-detecting it can be dangerous and the user may want
        # to define their own

#        if "ACE" not in atomsel().get('resname'):
#            minid = min(atomsel().get('resid'))
#            minresname = atomsel('resid %d' % minid).get('resname')[0]
#            print("\n\nINFO: Didn't find a C-terminal ACE for segment beginning"
#                  " with %s%d" % (minresname, minid))
#            patches += self._get_patch(seg, minid)
#
#        maxid = max(atomsel().get('resid'))
#        maxresname = atomsel('resid %d' % maxid).get('resname')
#        if "NMA" not in maxresname:
#            print("RESNAMES ARE")
#            print(set(atomsel().get('resname')))
#            print(set(atomsel().get('resid')))
#            print("\n\nINFO: Didn't find a N-terminal NMA for segment ending"
#                  " with %s%d" % (maxresname, maxid))
#            patches += self._get_patch(seg, maxid)

        # Now protein
        filename = self.tmp_dir + '/psf_protein_%s.pdb'% seg
        self._write_ordered_pdb(filename, 'all', prot_molid)
        print("Wrote %d atoms to the protein segment %s"
              % (len(atomsel('all')), seg))

        molecule.set_top(old_top)

        # Now write to psfgen input file
        string = '''
        pdbalias residue CYX CYS
        set protnam %s
        segment %s {
          first none
          last none
          pdb $protnam
        } 
        ''' % (filename, seg)
        self.file.write(string)
        self.file.write(patches)
        self.file.write('   coordpdb $protnam %s\n' % seg)

        return filename

    #==========================================================================

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
        atomsel('water').set('resname', 'TIP3')
        atomsel('resname TIP3').set('chain', 'W')
        atomsel('resname TIP3 and element O').set('name', 'OH2')

        # Dowser can name water hydrogens strangely
        atomsel('resname TIP3 and name HW1').set('name', 'H1')
        atomsel('resname TIP3 and name HW2').set('name', 'H2')

        # Select all the waters. We'll use the user field to track which
        # ones have been written
        allw = atomsel('water and user 1.0')
        print("Found %d water residues" % len(set(atomsel('water and user 1.0').get('residue'))))
        num_written = len(allw)/(9999*3)+1
        print("Going to write %d files for %d water atoms"
              % (num_written, len(allw)))

        # Pull out and write 10k waters at a time
        for i in range(num_written):
            residues = list(set(allw.get('residue')))[:9999]
            batchtxt = 'residue ' + \
                       ' '.join([str(s) for s in set(residues)])
            batch = atomsel(batchtxt)
            try:
                batch.set('resid', [k for k in range(1, len(batch)/3+1)
                                    for _ in range(3)])
            except ValueError:
                print("\nERROR! You have some waters missing hydrogens!\n"
                      "Found %d water residues, but %d water atoms. Check "
                      " your crystallographic waters in the input structure."
                      % (len(residues), len(batch)))
                quit(1)
            temp = tempfile.mkstemp(suffix='_%d.pdb' % i, prefix='psf_wat_',
                                    dir=self.tmp_dir)[1]
            batch.set('user', 0.0)
            batch.write('pdb', temp)
            allw.update()

        string = '''
        set waterfiles [glob -directory %s psf_wat_*.pdb]
        set i 0
        foreach watnam $waterfiles {
           segment W${i} {
              auto none
              first none
              last none
              pdb $watnam
           }
           coordpdb $watnam W${i}
           incr i
        }
          ''' % self.tmp_dir
        self.file.write(string)
        molecule.set_top(old_top)

        return num_written

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
        residues = list(set(alll.get('residue')))
        residues.sort()

        # Sanity check for < 10k lipids
        if len(residues) >= 10000:
            raise NotImplementedError("More than 10k lipids found")

        # Loop through all residues and renumber and correctly name them
        counter = 1
        for res in residues:
            # Renumber residue
            atomsel('residue %d' % res).set('resid', counter)
            counter = counter + 1 

            # Pick which naming dictionary to use
            resname = set(atomsel('residue %d' % res).get('resname'))
            if len(resname) != 1:
                raise ValueError("More than one name for residue %d" % res)
            resname = resname.pop()
            names = {}
            if resname == "POPC":
                # Carbons and hydrogen above nitrogen in head group
                # These selections must be done in order because there is a swap
                atomsel('(%s) and name C11' % self.lipid_sel).set('name', 't12')
                atomsel('(%s) and name C15' % self.lipid_sel).set('name', 'C11')
                atomsel('(%s) and name C14' % self.lipid_sel).set('name', 'C15')
                atomsel('(%s) and name C12' % self.lipid_sel).set('name', 'C14')
                atomsel('(%s) and name t12' % self.lipid_sel).set('name', 'C12')

                popc_names = {'H31':'H13A', 'H32':'H13B', 'H33':'H13C',
                              'H41':'H15A', 'H42':'H15B', 'H43':'H15C',
                              'H21':'H14A', 'H22':'H14B', 'H23':'H14C',
                              'H51':'H11A', 'H52':'H11B',
                              'H11':'H12A', 'H12':'H12B',
                              # Phosphate and its oxygens
                              'P1' :  'P', 'O1' :'O12', 'O2' :'O11',
                              'O3' :'O13', 'O4' :'O14'}
                names = popc_names
            elif resname == "POPE":
                pass
            elif resname == "POPG":
                pass
            else:
                raise NotImplementedError("Lipid %s unsupported" % resname)

            for name in names:
                atomsel('residue %d and name %s' 
                        % (res, name)).set('name', names[name])

        # Write temporary lipid pdb
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_lipid_',
                                dir=self.tmp_dir)[1]
        alll.set('user', 0.0)
        alll.write('pdb', temp)

        # Write to file
        string = '''
       set lipidfile %s
       set mid [mol new $lipidfile]
       segment L {
          first none
          last none
          pdb $lipidfile
       }
       coordpdb $lipidfile L
       mol delete $mid
        ''' % temp
        self.file.write(string)

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
        not_ions = atomsel('(same fragment as element Na Cl K)' \
                           ' and (not index %s)' % ' '.join([str(s) for s in set(total.get('index'))]))
        ions = set(total.get('residue')) - set(not_ions.get('residue'))

        if not len(ions): return
        ionstr = "residue " + " ".join([str(s) for s in ions])

        # Fix the names
        atomsel('%s and name NA' % ionstr).set('name', 'SOD')
        atomsel('%s and name CL' % ionstr).set('name', 'CLA')
        atomsel('%s and name K' % ionstr).set('name', 'POT')
        atomsel('%s and name NA' % ionstr).set('resname', 'SOD')
        atomsel('%s and name CL' % ionstr).set('resname', 'CLA')
        atomsel('%s and name K' % ionstr).set('resname', 'POT')

        # Renumber the residues since some may be above 10k
        residues = atomsel('name SOD CLA POT').get('residue')
        batch = atomsel('residue %s' % ' '.join([str(s) for s in set(residues)]))
        batch.set('resid', [k for k in range(1, len(batch)+1)])

        # Save the temporary ions file
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_ions_',
                                dir=self.tmp_dir)[1]
        atomsel('name SOD CLA POT').set('user', 0.0) # mark as saved
        atomsel('name SOD CLA POT').write('pdb', temp)

        string = '''
       set ionfile %s
       segment I {
          pdb $ionfile 
          first none
          last none
       }
       coordpdb $ionfile I
        ''' % temp
        self.file.write(string)
        molecule.set_top(old_top)

    #==========================================================================

    def _write_ligand_blocks(self, residue):
        """
        Matches ligands to available topology file, renames atoms, and then
        writes temporary files for the ligands

        Args:
          residue (int): Residue number of the ligand that will be written

        Returns:
          True if successful
        """
        # Put our molecule on top to simplify atom selection language
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        alig = atomsel('user 1.0 and residue %d' % residue)
        # Get a residue name charmm knows about, either through
        # manual translation or prompting the user
        res = alig.get('resname')[0]
        while not self._find_residue_in_rtf(residue=residue, molid=self.molid):
            print("\nERROR: Residue name %s wasn't found in any input "
                  "topology.\nWould you like to rename it?\n" % res)
            sys.stdout.flush()
            newname = raw_input("New residue name or CTRL+D to quit > ")
            sys.stdout.flush()
            alig.set('resname', newname)
            res = newname

        # Write temporary file containg the ligand and update tcl commands
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_ligand_',
                                dir=self.tmp_dir)[1]
        string = '''
       set ligfile %s
       segment %s {
         pdb $ligfile 
         first none
         last none
       }
       coordpdb $ligfile %s
        ''' % (temp, alig.get('resname')[0], alig.get('resname')[0])
        alig.write('pdb', temp)
        alig.set('user', 0.0)
        self.file.write(string)
        molecule.set_top(old_top)

    #==========================================================================

    def _write_ordered_pdb(self, filename, sel, molid=0):
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
        # Much much faster then get resids
        resids = set(atomsel(sel).get('residue'))
        idx = 1
        # For renumbering capping groups
        for rid in resids:
            # Handle bug where capping groups in same residue as the
            # neighboring amino acid Maestro writes it this way for some
            # reason but it causes problems down the line when psfgen doesn't
            # understand the weird combined residue
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
                    if len(atomsel('(%s) and resid %d' % (sel, resid-1))):
                        raise ValueError('ACE resid collision number %d' % resid-1)
                    atomsel('residue %d and resname ACE'
                            % rid).set('resid', resid-1)

                    # Handle all of the ACE atoms before the others
                    atoms = atomsel('residue %d and resname ACE'
                                    % rid).get('index')
                    atoms.extend(atomsel('residue %d and not resname ACE'
                                         % rid).get('index'))

                elif 'NMA' in names:
                    # Set NMA residue number as one more
                    resid = atomsel('residue %d and not resname NMA' % rid).get('resid')[0]
                    if len(atomsel('(%s) and resid %d' % (sel, resid+1))):
                        raise ValueError('NMA resid collision number %d' % resid+1)

                    atomsel('residue %d and resname NMA'
                            % rid).set('resid', resid+1)
                    # Handle all the NMA atoms after the others
                    atoms = atomsel('residue %d and not resname NMA'
                                    % rid).get('index')
                    atoms.extend(atomsel('residue %d and resname NMA'
                                         % rid).get('index'))

            else:
                atoms = atomsel('residue %d' % rid).get('index')

            for i in atoms:
                a = atomsel('index %d' % i)

                entry = ('%-6s%5d %-5s%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f'
                         '     %-4s%2s\n' % ('ATOM', idx, a.get('name')[0],
                                             a.get('resname')[0],
                                             a.get('chain')[0],
                                             a.get('resid')[0], 
                                             a.get('x')[0],
                                             a.get('y')[0],
                                             a.get('z')[0],
                                             0.0, 0.0, a.get('segname')[0],
                                             a.get('element')[0]))
                idx += 1
                fileh.write(entry)
        fileh.write('END\n')
        atomsel(sel).set('user', 0.0) # Mark as written
        fileh.close()
        molecule.set_top(old_top)


    #==========================================================================

    def _check_psf_output(self):
        """
        Scans the output psf from psfgen for atoms where the coordinate
        could not be set, indicating an unmatched atom. This checek is necessary
        because sometimes psfgen will run with no errors or warnings but will
        have unmatched atoms that are all at (0,0,0).
        """

        # Check file was written at all
        if not os.path.isfile('%s.pdb'% self.psf_name):
            print("\nERROR: psf file failed to write.\n"
                  "       Please see log above.\n")
            quit(1)

        problem_lines = []
        fileh = open('%s.pdb'% self.psf_name, 'r')
        for line in iter(fileh):
            # Use rsplit since some fields near beginning can mush together
            # 3rd field from end not fourth since element name will be absent
            words = line.rsplit()
            if words[0] == 'ATOM' and words[-3] == '-1.00':
                problem_lines.append(line[:-1])
        fileh.close()

        # Print out error messages
        if len(problem_lines):
            print("\nERROR: Couldn't find the following atoms.\n"
                  "Check if they are present in the original structure. "
                  "If they are, check dabble name translation or file a "
                  "bug report to Robin.\n")
            print('\n'.join(problem_lines))
            quit(1)
        else:
            print("\nINFO: Checked output pdb/psf has all atoms present "
                  "and correct.\n")


    #==========================================================================

    def _get_atoms_from_rtf(self, text, resname):
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
        for i in range(len(text)):
            words = text[i].split()
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

    def _find_residue_in_rtf(self, residue, molid):
        """
        Scans the input topology files to find a name match for the given
        residue ID, then pulls out the atoms involved and checks that they
        are all present in the input coordinates, prompting the user to correct
        the names of atoms that could not be matched.

        Residue ID is used because there can be multiple copies of a residue
        with the same name, but only one has missing or extra atoms.

        Args:
          topologies (list of str): Filenames of topologies to search
          residue (int): Residue number to search for (not resid)
          molid (int): VMD molecule ID

        Returns:
          True if all matching was successful
          False if the residue name cannot be found
        """

        resname = atomsel('residue %s and user 1.0'
                          % residue, molid=molid).get('resname')[0]
        print("Finding residue name %s" % resname)
        for top in self.topologies:
            topfile = open(top, 'r')
            topo_atoms = self._get_atoms_from_rtf(text=topfile.readlines(),
                                                  resname=resname)
            # Use first definition found of this residue
            if len(topo_atoms):
                break
            topfile.close()
        if not len(topo_atoms):
            return False
        print("INFO: Successfully found residue %s in input topologies" % resname)

        # Match up atoms with python sets
        pdb_atoms = set(atomsel('residue %s and user 1.0'
                                % residue, molid=molid).get('name'))
        pdb_only = pdb_atoms - topo_atoms
        topo_only = topo_atoms - pdb_atoms

        # If uneven number of atoms, there are missing or additional atoms
        if len(pdb_atoms) > len(topo_atoms):
            print("\nERROR: Cannot process modified residue %s.\n"
                  "       There are %d extra atoms in the input structure "
                  "that are undefined in the topology file. The "
                  "following atoms could not be matched and may "
                  "either be misnamed, or additional atoms. Please "
                  "check your input."
                  % (resname, len(pdb_atoms)-len(topo_atoms)))
            print("       [ %s ]\n" % ' '.join(pdb_only))
            print("       Cannot continue.\n")
            quit(1)
        if len(topo_atoms) > len(pdb_atoms):
            print("\nERROR: Cannot process modified residue %s.\n"
                  "       There are %d missing atoms in the input structure "
                  " that are defined in the topology file. The "
                  " following atoms could not be matched and may "
                  " either be misnamed or deleted atoms. Please "
                  " check your input."
                  % (resname, len(topo_atoms)-len(pdb_atoms)))
            print("       [ %s ]\n" % ' '.join(topo_only))
            print("       Cannot continue.\n")
            print("Found is %s\n" % pdb_atoms)
            quit(1)

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

                newname = raw_input("  %s  -> " % unmatched)
                while newname not in topo_only:
                    print("'%s' is not an available name in the topology."
                          "Please try again.\n" % newname)
                    newname = raw_input("  %s  -> " % unmatched)
                atomsel('residue %s and user 1.0 and name %s'
                        % (residue, unmatched)).set('name', newname)
                pdb_atoms = set(atomsel('residue %s and user 1.0'
                                        % residue).get('name'))
                topo_only = topo_atoms-pdb_atoms
                resname = newname

            # Recurse to check that everything is assigned correctly
            self._find_residue_in_rtf(residue, molid)
        print("INFO: Matched up all atom names for residue %s\n" % resname)
        return True

    #==========================================================================

    def _set_disulfide_bridges(self):
        """
        Adds PATCH lines corresponding to disulfide bridges.
        """

        # Simplify atom selection by putting relevant molecule on top
        old_top = molecule.get_top()
        molecule.set_top(self.molid)

        patches = ""
        indices = set(atomsel('name SG and resname CYX').get('index'))
        indices.update(atomsel('name SG and resname CYS '
                               'and not same residue as name HG').get('index'))
        while len(indices) > 0:
            idx1 = indices.pop()
            matches = set(atomsel('name SG and (resname CYX or (resname CYS '
                                  'and not same residue as name HG)) '
                                  'and not index %d and within 2.5 of index %d'
                                  % (idx1, idx1)))

            # Sanity check
            if len(matches) > 1:
                print('\nERROR: Found more than one possible disulfide bond '
                      'partner for atom %d. Don\'t know what to do now... '
                      'quack quack quack goodbye' % idx1)
                quit(1)
            elif len(matches) < 1:
                print('\nERROR: Found no disulfide bond partner for atom %d '
                      'Please check the input file is prepared properly' % idx1)
                quit(1)
            idx2 = matches.pop()
            indices.remove(idx2)

            # Set up the disu patch line
            res1 = atomsel('index %d' % idx1).get('resid')[0]
            res2 = atomsel('index %d' % idx2).get('resid')[0]
            seg1 = atomsel('index %d' % idx1).get('fragment')[0]
            seg2 = atomsel('index %d' % idx2).get('fragment')[0]

            print("\nINFO: Disulfide bond between residues %d and %d"
                  % (res1, res2))
            patches += 'patch DISU P%s:%d P%s:%d\n' % (seg1, res1, seg2, res2)

        self.file.write(patches)
        molecule.set_top(old_top)

    #========================================================================== 

    def _get_bonded_atoms(self, molid, index):
        """
        Returns the element of all atoms bonded to the current atom.

        Args:
           molid (int): VMD molecule ID to consider
           index (int): Atom index to look at bonded atoms

        Returns:
          (list of str) elements of atoms bound to the current atom
        """

        a = atomsel('index %d' % index, molid=molid)
        bound = []
        for atom in a.bonds[0]:
            bound.append(atomsel('index %d' % atom).get('element')[0])
        return bound
    
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
        print("Or type HELP for a list of all patches I know about")
        patchname = raw_input("> ")
        if patchname == "HELP":
            print("   PATCH     COMMENT")
            print("   -----     -------")
            for p in avail_patches:
                print("%7s %s" % (p, avail_patches[p]))
            patchname = raw_input("> ")
        while patchname not in avail_patches:
            print("I don't know about patch %s" % patchname)
            patchname = raw_input("Try again > ")

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
                if not len(tokens): continue 
                if tokens[0] == "PRES": 
                    comment = ' '.join(tokens[tokens.index("!")+1:])
                    avail_patches[tokens[1]] = comment
        return avail_patches

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
