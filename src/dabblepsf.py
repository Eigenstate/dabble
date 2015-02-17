import molecule
import glob,os,tempfile
from atomsel import *

# Writes a psf file with the given options
def write_psf(psf_name, molid=0, lipid_sel="lipid"):
    # Clean up all temp files from previous runs
    # TODO: Move to end if successful

    print("Removing %d old psfgen blocks" % len(glob.glob("/tmp/psf_*.pdb")))
    for oldfile in glob.glob("/tmp/psf_*.pdb") :
        os.remove(oldfile)
    for oldfile in glob.glob("/tmp/dabble_psfgen*.tcl") :
        os.remove(oldfile)

    filename = tempfile.mkstemp(suffix='.tcl', prefix='dabble_psfgen')[1]
    file = open(filename,'w')

    #  Finds the psfgen package and sets the output file name
    string='''
    set dir [file join $env(VMDDIR) plugins [vmdinfo arch] tcl psfgen1.6]
    package ifneeded psfgen 1.6.2 [list load [file join $dir libpsfgen.so]]
    package require psfgen
    set output "%s"
    resetpsf
    set dabbledir $::env(DABBLEDIR)
    ''' % psf_name
    file.write(string)

    # Specify default parameter sets
    string='''
    topology ${dabbledir}/charmm_params/top_water_ions.rtf
    topology ${dabbledir}/charmm_params/top_all36_cgenff.rtf
    topology ${dabbledir}/charmm_params/top_all36_prot.rtf
    topology ${dabbledir}/charmm_params/top_all36_lipid.rtf
    '''
    file.write(string)

    # Ask the user what parameter sets they want to use
    # TODO: Only ask if a ligand or something is found
    print("""\nWhich CHARMM topology files should I use?"
    Currently using:
       - top_water_ions.rtf
       - top_all36_cgenff.rtf
       - top_all36_prot.rtf
       - top_all36_lipid.rtf """)
    print("Enter the full path to the filename(s), separated by a comma, of any "
          "other rtf files you wish to use.\n")
    inp = raw_input('> ')
    if inp:
        rtfs = inp.split(',')
        for t in rtfs: 
            file.write("topology %s\n" % t)

    # Save water 10k molecules at a time
    if len(atomsel('water',molid=molid)) > 0:
        write_water_blocks(molid=molid)
        string = '''
        set waterfiles [glob -directory /tmp psf_wat_*.pdb]
        set i 0
        foreach watnam $waterfiles {
          set mid [mol new $watnam]
          set water [atomselect top "water"]
          segment W${i} {
            auto none
            pdb $watnam
          }
          pdbalias atom HOH O OH2
          coordpdb $watnam W${i}
          mol delete $mid
          incr i
        }
        '''
        file.write(string)
    # End water

    # Now ions if present, changing the atom names
    if len(atomsel('ions',molid=molid)) > 0 :
        temp = write_ion_blocks(molid=molid)
        string='''
        set ionfile %s
        set mid [mol new $ionfile]
        segment I {
          pdb $ionfile 
        }
        coordpdb $ionfile I
        mol delete $mid ''' % temp
        file.write(string)
    # End ions

    # Now lipid
    if len(atomsel(lipid_sel)) > 0:
        temp=write_lipid_blocks(lipid_sel=lipid_sel,molid=molid)
        string='''
          set lipidfile %s
          set mid [mol new $lipidfile]
          segment L {
            pdb $lipidfile
          }
          coordpdb $lipidfile L
          mol delete $mid
        ''' % temp
        file.write(string)
    # End lipid

    # Save the protein with correct atom names
    # TODO: Multiple chains
    if len(atomsel('protein',molid=molid)) > 0:
        temp = write_protein_blocks(file,molid=molid)
        # End protein

    # Check if there is anything else and let the user know about it
    leftovers = atomsel('not water and not lipid and not (protein or resname ACE) and not ions',molid=molid)
    if len(leftovers) > 0:
        temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_extra_')[1]
        leftovers.write('pdb',temp)
        print("\n\nALERT: Found what is probably a ligand!")
        print("Residue names are: ")
        print set(leftovers.get('resname'))
        print("Saving to '%s'" % temp)
        print("Good luck with psfgen lololol")

    #protein_opts='\n'.join(iter(raw_input, 'END'))

    # Write the output files and run
    string='''
    writepsf x-plor cmap ${output}.psf
    writepdb ${output}.pdb'''
    file.write(string)
    file.close()

    from VMD import evaltcl
    evaltcl('play %s' % filename)

# Writes the temporary protein PDB file with correct atom names for psfgen
# This is also the worst.
def write_protein_blocks(file,molid=0):
    # Put molid on top to simplify atom selections
    old_top=molecule.get_top()
    molecule.set_top(molid)

    # Check the protein has hydrogens
    if len(atomsel('protein and element H')) == 0:
        print("\n\nERROR: There are no hydrogens on your protein")
        print("       You need to add them before parameterizing")
        print("       because psfgen cannot be trusted to do it correctly.\n")
        quit(1)

    # Terminal residues ACE TODO NMA
    atomsel('resname ACE and name CH3').set('name','CAY')
    atomsel('resname ACE and name HH31').set('name','HY1')
    atomsel('resname ACE and name HH32').set('name','HY2')
    atomsel('resname ACE and name HH33').set('name','HY3')
    atomsel('resname ACE and name C').set('name','CY')
    atomsel('resname ACE and name O').set('name','OY')

    # Disulfide briges
    disulfide_patches = ''
    indices = set( atomsel('name SG and resname CYX').get('index') )
    while len(indices) > 0 :
        idx1 = indices.pop()
        matches = set( atomsel('name SG and resname CYX and not index %d and within 2.5 of index %d' % (idx1, idx1) ))

        # Sanity check
        if len(matches) > 1 :
            print('\nERROR: Found more than one possible disulfide bond partner for atom %d'
                    '       Don\'t know what to do now... quack quack quack goodbye' % idx1)
            quit(1)
        elif len(matches) < 1 :
            print('\nERROR: Found no disulfide bond partner for atom %d'
                    '       Please check your input file is prepared properly' % idx1)
        idx2 = matches.pop()
        indices.remove(idx2)

        # Set up the disu patch line
        res1 = atomsel('index %d' % idx1).get('resid')[0]
        res2 = atomsel('index %d' % idx2).get('resid')[0]
        print("\nINFO: Disulfide bond between residues %d and %d" % (res1,res2))
        disulfide_patches += 'patch DISU P:%d P:%d\n' % (res1,res2)

    # Rename CYX -> CYS for naming conventions
    atomsel('resname CYX').set('resname','CYS')

    # Histidine naming convention
    atomsel('resname HID').set('resname','HSD')
    atomsel('resname HIE').set('resname','HSE')
    atomsel('resname HIP').set('resname','HSP') 

    # Determine protonation states of those just called HIS based on atom names
    atomsel('resname HIS and same residue as name HE2').set('resname','HSE')
    atomsel('resname HIS and same residue as name HD1').set('resname','HSD')
    # NOTE: This atomsel MUST come after the previous two!
    atomsel('resname HIS and same residue as name HE2 and same residue as name HD1').set('resname',
            'HSP')
    
    # Isoleucine
    atomsel('resname ILE and name CD1').set('name','CD')
    atomsel('resname ILE and name HD11').set('name','HD1')
    atomsel('resname ILE and name HD12').set('name','HD2')
    atomsel('resname ILE and name HD13').set('name','HD3')
    atomsel('resname ILE and name HG12').set('name','HG11')
    atomsel('resname ILE and name HG13').set('name','HG12')

    # Glutamine check all residues for protonated one
    t=atomsel('resname GLU or resname GLH').get('residue')
    for residue in t:
        if "HE1" in atomsel('residue %s' % residue).get('name') :
            atomsel('residue %s and name HE1' % residue).set('name','HE2')
            atomsel('residue %s and name OE1' % residue).set('name','temp')
            atomsel('residue %s and name OE2' % residue).set('name','OE1')
            atomsel('residue %s and name temp' % residue).set('name','OE2')
            atomsel('residue %s').set('resname','GLUP')
        elif "HE2" in atomsel('residue %s' % residue).get('name') :
            atomsel('residue %s' % residue).set('resname','GLUP')
        # TODO: Is this correct? Atom name messups
        elif "HXT" in atomsel('residue %s' % residue).get('name') :
            atomsel('residue %s' % residue).set('resname','GLUP')
            atomsel('residue %s and name HXT' % residue).set('name','HE2')
        else :
            atomsel('residue %s' % residue).set('resname','GLU')

    # Aspartate check each residue to see if it's protonated
    t=atomsel('resname ASP or resname ASH').get('residue')
    for residue in t :
        if "HD1" in atomsel('residue %s' % residue).get('name') :
            atomsel('residue %s and name HD1' % residue).set('name','HD2')
            atomsel('residue %s and name OD1' % residue).set('name','temp')
            atomsel('residue %s and name OD2' % residue).set('name','OD1')
            atomsel('residue %s and name temp' % residue).set('name','OD2')
            atomsel('residue %s').set('resname','ASPP')
        if "HD2" in atomsel('residue %s' % residue).get('name') :
            atomsel('residue %s' % residue).set('resname','ASPP')
        else :
            atomsel('residue %s' % residue).set('resname','ASP')

    # Hydrogens
    atomsel('(resname ACE or protein) and name H').set('name','HN')
    atomsel('(resname ACE or protein) and name HA2').set('name','HA1')
    atomsel('(resname ACE or protein) and name HA3').set('name','HA2')
    atomsel('(resname ACE or protein) and name HG2').set('name','HG1')
    atomsel('(resname ACE or protein) and name HG3').set('name','HG2')
    atomsel('(resname ACE or protein) and name HB2 and not resname ALA').set('name','HB1')
    atomsel('(resname ACE or protein) and name HB3 and not resname ALA').set('name','HB2')
    atomsel('(resname ACE or protein) and name HD2 and not (resname HSP HSE HSD ASPP PHE)').set('name','HD1')
    atomsel('(resname ACE or protein) and name HD3 and not (resname HIS HSE HSD PHE)').set('name','HD2')
    atomsel('(resname ACE or protein) and name HE2 and not (resname TRP HSP HSE HSD GLUP PHE)').set('name','HE1')
    atomsel('(resname ACE or protein) and name HE3 and not (resname TRP HSP HSE HSD PHE)').set('name','HE2')
    atomsel('(resname ACE or protein) and name HG and resname SER CYS').set('name','HG1')

    # Renumber residues
    T = atomsel('resname ACE or protein')
    T.set('resid',T.get('residue'))

    # Now protein
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_protein_')[1]
    atomsel('resname ACE or protein').write('pdb',temp)
    print("Wrote %d atoms to the protein temp file" % len(atomsel('resname ACE or protein')))
    molecule.set_top(old_top)

    # Now write to psfgen input file
    string='''
          set protnam %s
          segment P {
            pdb $protnam
          } 
          ''' % temp
    file.write(string)
    file.write(disulfide_patches)
    file.write('coordpdb $protnam P') 

    return temp

# Writes a bunch of temp files with 10000 waters each, to bypass psfgen
# being stupid with files with more than 10000 of a residue
def write_water_blocks(molid=0):
    # Select all the waters. We'll use the user field to track which ones have been written
    atomsel('water',molid=molid).set('user', 0.0)
    all=atomsel('water and user 0.0',molid=molid)
    num_written = len(all)/(9999*3)+1
    print("Going to write %d files for %d atoms" % (num_written, len(all)))

    # Pull out and write 10k waters at a time
    for i in range(num_written) :
        residues = all.get('residue')[:9999*3] # 9999 molecules, 3 atoms each
        batch = atomsel('residue ' + ' '.join(map(str,set(residues))), molid=molid)
        try :
            batch.set('resid', [k for k in range(1,len(batch)/3+1) for _ in range(3)])
        except ValueError:
            print("ERROR! You have some waters missing hydrogens!")
            print("       Found %d water residues, but %d water atoms" %
                (len(residues),len(batch)))
            print("       Check your crystallographic waters in the input structure.\n")
            quit(1)
        temp = tempfile.mkstemp(suffix='_%d.pdb' % i, prefix='psf_wat_')[1]
        batch.write('pdb', temp)
        batch.set('user', 1.0)
        all.update()
    return num_written

# Gotta renumber the lipid residues too because some can have **** instead of resid
def write_lipid_blocks(lipid_sel="lipid",molid=0):
    # Set the user field of all lipids, we'll use it to keep track of which
    # were written
    atomsel(lipid_sel,molid=molid).set('user',0.0)
    all = atomsel('(%s) and user 0.0' % lipid_sel)

    # Detect what lipid kind this is. HACK for now, it really should figure it out,
    # but this is what works
    lipid_types = set( all.get('resname') )
    atoms_per_lip = 0
    lip_type = lipid_types.pop() # Should remove only element in this list

    if len(lipid_types) > 0:
      print("\nERROR: Found more than one type of lipid in your membrane.\n"
              "       This currently isn't supported in psfgen")
      quit(1)
    elif lip_type == 'POPC' :
      atoms_per_lip = 134
    else :
      print("\nERROR: Unsupported lipid type %s" % lipid_types[0])
      print("       Get Robin to put it in the code, it's quick.")
      quit(1)

    # Sanity check for correct number of atoms
    if len(all) % atoms_per_lip != 0:
      print("\nERROR: Incorrect number of atoms for these lipids\n"
              "       Have %d lipid atoms / %d atoms per lipid, "
              "       which results in %f lipids!" % (len(all), atoms_per_lip,
              len(all)/atoms_per_lip))
      quit(1)

    # Sanity check for <10k lipids
    num_written = len(all)/(atoms_per_lip)+1
    if num_written >= 10000 :
      print("\nERROR: Have more than 10k lipids\n"
              "       Support not implemented. Ask Robin to fix it.")
      quit(1)

    # Renumber resids starting from 0
    for i in range(num_written) :
        residues = all.get('residue')
        batch = atomsel('residue ' + ' '.join(map(str,set(residues))), molid=molid)
        try :
            batch.set('resid', [k for k in range(1,len(batch)/atoms_per_lip+1) for _ in range(atoms_per_lip)])
        except ValueError:
            print("\nERROR! Something went wrong renumbering the lipids!")
            quit(1)

    # Rename atoms. Wow, this is really the worst
    # Carbons and hydrogen above nitrogen in head group
    atomsel('(%s) and name C11' % lipid_sel).set('name','t12')
    atomsel('(%s) and name C15' % lipid_sel).set('name','C11')
    atomsel('(%s) and name C14' % lipid_sel).set('name','C15')
    atomsel('(%s) and name C12' % lipid_sel).set('name','C14')
    atomsel('(%s) and name t12' % lipid_sel).set('name','C12')

    atomsel('(%s) and name H31' % lipid_sel).set('name','H13A')
    atomsel('(%s) and name H32' % lipid_sel).set('name','H13B')
    atomsel('(%s) and name H33' % lipid_sel).set('name','H13C')
    atomsel('(%s) and name H41' % lipid_sel).set('name','H15A')
    atomsel('(%s) and name H42' % lipid_sel).set('name','H15B')
    atomsel('(%s) and name H43' % lipid_sel).set('name','H15C')
    atomsel('(%s) and name H21' % lipid_sel).set('name','H14A')
    atomsel('(%s) and name H22' % lipid_sel).set('name','H14B')
    atomsel('(%s) and name H23' % lipid_sel).set('name','H14C')

    atomsel('(%s) and name H51' % lipid_sel).set('name','H11A')
    atomsel('(%s) and name H52' % lipid_sel).set('name','H11B')
    atomsel('(%s) and name H11' % lipid_sel).set('name','H12A')
    atomsel('(%s) and name H12' % lipid_sel).set('name','H12B')

    # Phosphate and its oxygens
    atomsel('(%s) and name P1' % lipid_sel).set('name','P')
    atomsel('(%s) and name O1' % lipid_sel).set('name','O12')
    atomsel('(%s) and name O2' % lipid_sel).set('name','O11')
    atomsel('(%s) and name O3' % lipid_sel).set('name','O13')
    atomsel('(%s) and name O4' % lipid_sel).set('name','O14')

    # Fortunately all the other atom names are correct.
    # TODO: This will probably be different for non-POPC lipids...

    # Write temporary lipid pdb
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_lipid_')[1]
    atomsel(lipid_sel,molid=molid).write('pdb', temp)
    return temp

def write_ion_blocks(molid=0) :
    # Fix the names
    atomsel('ions and name NA').set('name','SOD')
    atomsel('ions and name CL').set('name','CLA')
    atomsel('ions and name K').set('name','POT')
    atomsel('ions and name NA').set('resname','SOD')
    atomsel('ions and name CL').set('resname','CLA')
    atomsel('ions and name K').set('resname','POT')

    # Renumber the residues since some may be above 10k
    residues = atomsel('ions',molid=molid).get('residue')
    batch = atomsel('residue ' + ' '.join(map(str,set(residues))), molid=molid)
    batch.set('resid', [k for k in range(1,len(batch)+1)])

    # Save the temporary ions file
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_ions_')[1]
    atomsel('ions',molid=molid).write('pdb',temp)
    return temp

