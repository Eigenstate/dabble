import molecule
import sys,glob,os,tempfile
from atomsel import *

# Writes a psf file with the given options
def write_psf(psf_name, molid=0, lipid_sel="lipid"):
    # Clean up all temp files from previous runs
    # TODO: Move to end if successful

    #print("Removing %d old psfgen blocks" % len(glob.glob("/tmp/psf_*.pdb")))
    #for oldfile in glob.glob("/tmp/psf_*.pdb") :
    #    os.remove(oldfile)
    #for oldfile in glob.glob("/tmp/dabble_psfgen*.tcl") :
    #    os.remove(oldfile)
    
    tmp_dir = tempfile.mkdtemp(prefix='dabble_psfgen')
    filename = tempfile.mkstemp(suffix='.tcl', prefix='dabble_psfgen', dir=tmp_dir)[1]
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
    topology ${dabbledir}/charmm_params/top_all36_caps.rtf
    topology ${dabbledir}/charmm_params/top_all36_cgenff.rtf
    topology ${dabbledir}/charmm_params/top_water_ions.rtf
    topology ${dabbledir}/charmm_params/top_all36_prot.rtf
    topology ${dabbledir}/charmm_params/top_all36_lipid.rtf
    topology ${dabbledir}/charmm_params/top_all36_carb.rtf
    '''
    file.write(string)

    # Ask the user what parameter sets they want to use
    # TODO: Only ask if a ligand or something is found
    sys.stdout.flush()
    print("""\nWhich CHARMM topology files should I use?"
    Currently using:
       - top_all36_caps.rtf (capping groups by Robin)
       - top_water_ions.rtf
       - top_all36_cgenff.rtf
       - top_all36_prot.rtf
       - top_all36_lipid.rtf
       - top_all36_carb.rtf""")
    print("Enter the full path to the filename(s), separated by a comma, of any "
          "other rtf files you wish to use.\n")
    sys.stdout.flush()
    inp = raw_input('> ')
    if inp:
        rtfs = inp.split(',')
        for t in rtfs: 
            file.write("topology %s\n" % t)

    # Mark all atoms as unsaved with the user field
    atomsel('all',molid=molid).set('user',1.0)

    # Now ions if present, changing the atom names
    if len(atomsel('element Na Cl K',molid=molid)) > 0 :
        temp = write_ion_blocks(file,tmp_dir,molid=molid)
    # End ions

    # Save water 10k molecules at a time
    if len(atomsel('water',molid=molid)) > 0:
        write_water_blocks(file,tmp_dir,molid=molid)
    # End water

    # Now lipid
    if len(atomsel(lipid_sel)) > 0:
        temp=write_lipid_blocks(file,tmp_dir,lipid_sel=lipid_sel,molid=molid)
    # End lipid

    # Here we redefine the protein segments, as if there are breaks in it with
    # capping groups at each end it can still be listed as the same segment by Maestro
    # but since it is not strictly connected should not be the same segment as psfgen 
    # will connect together all atoms in a linear chain.
    cap_n = sorted(set(atomsel('resname ACE').get('resid')))
    cap_c = sorted(set(atomsel('resname NMA').get('resid')))
    assert len(cap_c)==len(cap_n), "Uneven number of capping groups!"

    # TODO: Sanity check that internal caps are next to each other in sequence
    segnum=1
    while len(cap_c) > 0 :
        nid = cap_n.pop(0)
        cid = cap_c.pop(0)
        assert nid < cid, "N-terminal residue before C-terminal? NMA resid %d, ACE resid %d"% (nid,cid)
        segment = atomsel('resid >= %d and resid <= %d'% (nid,cid))
        segment.set('segname','P%s' % segnum)
        segment.set('segid','P%s' % segnum)
        # This is stupid, but to fix the Maestro issue with combined capping group and regualr
        # residues, we have to renumber the residues, save, and reload so that VMD splits the 
        # combined residue into 2
        print("Writing protein chain from residue %d to %d"% (nid,cid))
        temp = tempfile.mkstemp(suffix='_P%s.pdb' % segnum, prefix='psf_prot_', dir=tmp_dir)[1]
        write_ordered_pdb(temp, sel='segname P%s'% segnum, molid=molid)

        prot_molid = molecule.load('pdb', temp)
        write_protein_blocks(file, seg='P%s'% segnum, tmp_dir=tmp_dir,molid=prot_molid)
        molecule.delete(prot_molid)
        #os.remove(temp)
        segnum += 1
    # End protein

    # Check if there is anything else and let the user know about it
    leftovers = atomsel('user 1.0', molid=molid)
    if len(leftovers) > 0 :
        print("Found extra ligands: %s" % set(leftovers.get('resname')))
    for l in set( leftovers.get('resid') ) :
        write_ligand_blocks(file, tmp_dir, resid=l, molid=molid)

    # Write the output files and run
    string='''
    writepsf x-plor cmap ${output}.psf
    writepdb ${output}.pdb'''
    file.write(string)
    file.close()

    from VMD import evaltcl
    evaltcl('play %s' % filename)
    check_psf_output(psf_name)

# Writes the temporary protein PDB file with correct atom names for psfgen
# This is also the worst.
def write_protein_blocks(file,tmp_dir,seg,molid=0):
    # Put molid on top to simplify atom selections
    old_top=molecule.get_top()
    molecule.set_top(molid)

    # Renumber residues starting from 1 TODO: This doesnt start from 1
    T = atomsel('segname %s' % seg)
    T.set('resid',[r+1 for r in T.get('residue')])

    # Check the protein has hydrogens
    if len(atomsel('segname %s and element H'% seg)) == 0:
        print("\n\nERROR: There are no hydrogens on your protein")
        print("       You need to add them before parameterizing")
        print("       because psfgen cannot be trusted to do it correctly.\n")
        quit(1)

    # Terminal residue ACE
    ace_names = {'CH3' :'CAY', 'O'   :'OY', 
                 'HH31':'HY1', 'HH32':'HY2', 'HH33':'HY3',
                 'H1'  :'HY1', 'H2'  :'HY2', 'H3'  :'HY3',
                 '1H'  :'HY1', '2H'  :'HY2', '3H'  :'HY3' }
    for n in ace_names :
        atomsel('segname %s and resname ACE and name %s'%(seg,n)).set('name',ace_names[n])

    # Terminal residue NMA
    nma_names = {'HN'  :'HNT', 'H'   :'HNT', 
                 'CT3' :'CAT', 'CA'  :'CAT', 'CH3' :'CAT',
                 'HA1' :'HT1', 'HA2' :'HT2', 'HA3' :'HT3',
                 '1HA' :'HT1', '2HA' :'HT2', '3HA' :'HT3', 
                 'HH31':'HT1', 'HH32':'HT2', 'HH33':'HT3' }
    for n in nma_names :
        atomsel('segname %s and resname NMA and name %s' %(seg,n)).set('name',nma_names[n])

    patches = ''

    # Disulfide briges
    indices = set( atomsel('segname %s and name SG and resname CYX'% seg).get('index') )
    while len(indices) > 0 :
        idx1 = indices.pop()
        matches = set( atomsel('segname %s and name SG and resname CYX and not index %d and within 2.5 of index %d' % (seg, idx1, idx1) ))

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
        patches += 'patch DISU %s:%d %s:%d\n' % (seg,res1,seg,res2)

    # Rename CYX -> CYS for naming conventions
    atomsel('resname CYX').set('resname','CYS')
    # Histidine naming convention
    atomsel('segname %s and resname HID'% seg).set('resname','HSD')
    atomsel('segname %s and resname HIE'% seg).set('resname','HSE')
    atomsel('segname %s and resname HIP'% seg).set('resname','HSP') 

    # Determine protonation states of those just called HIS based on atom names
    atomsel('segname %s and resname HIS and same residue as name HE2'% seg).set('resname','HSE')
    atomsel('segname %s and resname HIS and same residue as name HD1'% seg).set('resname','HSD')
    # NOTE: This atomsel MUST come after the previous two!
    atomsel('segname %s and resname HIS and same residue as name HE2 and same residue as name HD1' % seg).set('resname','HSP')
    
    # Isoleucine
    iso_names = { 'CD1' :  'CD',
                  'HD11': 'HD1', 'HD12':'HD2', 'HD13':'HD3',
                  'HG12':'HG11', 'HG13':'HG12' }
    for n in iso_names :
        atomsel('segname %s and resname ILE and name %s' %(seg,n)).set('name',iso_names[n])
                  
## TODO DEBUG bug in here for naomi's protein
    # Glutamine check all residues for protonated one
    # Can't use dictionary for names since two of them must be swapped
    # Must loop by resid, not residue, to get correct PATCH statement index
# TODO: Is this correct? Atom name messups
    t = set( atomsel('segname %s and (resname GLU GLH GLUP)'% seg).get('resid') )
    for resid in t:
        atomsel('segname %s and resid %s' % (seg,resid)).set('resname','GLU')
        if "HE1" in atomsel('resid %s' % resid).get('name') :
            atomsel('segname %s and resid %s and name HE1'% (seg,resid)).set('name','HE2')
            atomsel('segname %s and resid %s and name OE1'% (seg,resid)).set('name','temp')
            atomsel('segname %s and resid %s and name OE2'% (seg,resid)).set('name','OE1')
            atomsel('segname %s and resid %s and name temp'% (seg,resid)).set('name','OE2')
            patches += 'patch GLUP %s:%d\n' % (seg, resid)
        elif "HE2" in atomsel('segname %s and resid %s' % (seg,resid)).get('name') :
            patches += 'patch GLUP %s:%d\n' % (seg, resid)
        elif "HXT" in atomsel('segname %s and resid %s' % (seg,resid)).get('name') :
            atomsel('segname %s and resid %s and name HXT' % (seg,resid)).set('name','HE2')
            patches += 'patch GLUP %s:%d\n' % (seg, resid)
## END BUG LOCATION

    # Aspartate check each residue to see if it's protonated
    # Can't use dictionary for names since two of them must be swapped
    # Must loop by resid, not residue, to get correct PATCH statement index
    t = set( atomsel('resname ASP ASH ASPP').get('resid') )
    for resid in t :
        atomsel('segname %s and resid %s'% (seg,resid)).set('resname','ASP')
        if "HD1" in atomsel('segname %s and resid %s' % (seg,resid)).get('name') :
            atomsel('segname %s and resid %s and name HD1' % (seg,resid)).set('name','HD2')
            atomsel('segname %s and resid %s and name OD1' % (seg,resid)).set('name','temp')
            atomsel('segname %s and resid %s and name OD2' % (seg,resid)).set('name','OD1')
            atomsel('segname %s and resid %s and name temp' % (seg,resid)).set('name','OD2')
            atomsel('segname %s and resid %s').set('resname','ASPP')
            patches += 'patch ASPP %s:%d\n' % (seg,resid)
        if "HD2" in atomsel('segname %s and resid %s' % (seg,resid)).get('name') :
            patches += 'patch ASPP %s:%d\n' % (seg,resid)

    # Hydrogens can be somewhat residue-dependent
    h_names = {'H'  : 'HN', 'HA2':'HA1',
               'HA3':'HA2', 'HG2':'HG1',
               'HG3':'HG2'}
    for h in h_names :
        atomsel('segname %s and name %s' % (seg,h)).set('name',h_names[h])

    # These two statements must execute in this specific order
    atomsel('segname %s and name HB2 and not (resname ALA)'% seg).set('name','HB1')
    atomsel('segname %s and name HB3 and not (resname ALA)'% seg).set('name','HB2')

    # Some residue-specific stuff
    atomsel('segname %s and name H2 and resname MET'% seg).set('name','HN')
    atomsel('segname %s and name HD2 and not (resname TYR ILE HSP HSE HSD ASP PHE)' % seg).set('name','HD1')
    atomsel('segname %s and name HD3 and not (resname TYR ILE HIS HSE HSD ASP PHE)' % seg).set('name','HD2')
    atomsel('segname %s and name HE2 and not (resname TYR MET TRP HSP HSE HSD GLUP PHE)' % seg).set('name','HE1')
    atomsel('segname %s and name HE3 and not (resname TYR MET TRP HSP HSE HSD PHE)' % seg).set('name','HE2')
    atomsel('segname %s and name HG and resname SER CYS' % seg).set('name','HG1')

    # Now protein
    filename = tmp_dir + '/psf_protein_%s.pdb'% seg
    write_ordered_pdb(filename, sel='segname %s' % seg, molid=molid)
    print("Wrote %d atoms to the protein segment %s" % (len(atomsel('segname %s' % seg)),seg))
    molecule.set_top(old_top)

    # Now write to psfgen input file
    string='''
    set protnam %s
    segment %s {
      first none
      last none
      pdb $protnam
    } 
    ''' % (filename,seg)
    file.write(string)
    file.write(patches)
    file.write('   coordpdb $protnam %s\n' % seg)
    #file.write('   coordpdb $protnam %s\nguesscoord' % seg) 

    return filename 

# Writes a bunch of temp files with 10000 waters each, to bypass psfgen
# being stupid with files with more than 10000 of a residue
def write_water_blocks(file,tmp_dir,molid=0):
    # Set consistent residue and atom names, crystal waters can be named HOH, etc
    atomsel('water not name CLA SOD POT').set('resname','TIP3')
    atomsel('water not name CLA SOD POT').set('chain','W')
    atomsel('chain W and name O').set('name','OH2')

    # Select all the waters. We'll use the user field to track which ones have been written
    all=atomsel('chain W and user 1.0',molid=molid)
    num_written = len(all)/(9999*3)+1
    print("Going to write %d files for %d water atoms" % (num_written, len(all)))

    # Pull out and write 10k waters at a time
    for i in range(num_written) :
        residues = list(set(all.get('residue')))[:9999] # 9999 unique residues
        batchtxt = 'residue ' + ' '.join(map(str,set(residues)))
        batch = atomsel(batchtxt, molid=molid)
        try :
            batch.set('resid', [k for k in range(1,len(batch)/3+1) for _ in range(3)])
        except ValueError:
            print("ERROR! You have some waters missing hydrogens!")
            print("       Found %d water residues, but %d water atoms" %
                (len(residues),len(batch)))
            print("       Check your crystallographic waters in the input structure.\n")
            quit(1)
        temp = tempfile.mkstemp(suffix='_%d.pdb' % i, prefix='psf_wat_', dir=tmp_dir)[1]
        batch.set('user', 0.0)
        batch.write('pdb',temp)
        all.update()

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
      ''' % tmp_dir
    file.write(string)

    return num_written

# Gotta renumber the lipid residues too because some can have **** instead of resid
def write_lipid_blocks(file,tmp_dir,lipid_sel="lipid",molid=0):
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
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_lipid_', dir=tmp_dir)[1]
    atomsel(lipid_sel,molid=molid).set('user', 0.0)
    atomsel(lipid_sel,molid=molid).write('pdb',temp)

    # Write to file
    string='''
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
    file.write(string)

    return temp


def write_ion_blocks(file,tmp_dir,molid=0) :
    # Fix the names
    atomsel('name NA').set('name','SOD')
    atomsel('name CL').set('name','CLA')
    atomsel('name K').set('name','POT')
    atomsel('name NA').set('resname','SOD')
    atomsel('name CL').set('resname','CLA')
    atomsel('name K').set('resname','POT')

    # Renumber the residues since some may be above 10k
    residues = atomsel('name SOD CLA POT',molid=molid).get('residue')
    batch = atomsel('residue ' + ' '.join(map(str,set(residues))), molid=molid)
    batch.set('resid', [k for k in range(1,len(batch)+1)])

    # Save the temporary ions file
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_ions_', dir=tmp_dir)[1]
    atomsel('name SOD CLA POT',molid=molid).set('user',0.0) # mark as saved
    atomsel('name SOD CLA POT',molid=molid).write('pdb',temp)

    string='''
   set ionfile %s
   segment I {
      pdb $ionfile 
      first none
      last none
   }
   coordpdb $ionfile I
    ''' % temp
    file.write(string)

    return temp


def write_ligand_blocks(file, tmp_dir, resid, molid=0) :
    A = atomsel('user 1.0 and resid %d' % resid)

    # Adjust residue name if known ligand
    warning = True
    if 'CLR' in A.get('resname') :
        print("INFO: Found a cholesterol!")
        warning = False
        # Reassign names
        clol_names = {'H11' : 'H1A',  'H12' : 'H1B',
                      'H21' : 'H2A',  'H22' : 'H2B',
                      'O1'  : 'O3',   'HO1' : 'H3\'', 
                      'H41 ': 'H4A',  'H42' : 'H4B',
                      'H71' : 'H7A',  'H72' : 'H7B',
                      'H111': 'H11A', 'H112': 'H11B',
                      'H121': 'H12A', 'H122': 'H12B',
                      'H151': 'H15A', 'H152': 'H15B',
                      'H161': 'H16A', 'H162': 'H16B',
                      'H181': 'H18A', 'H182': 'H18B', 'H183': 'H18C',
                      'H191': 'H19A', 'H192': 'H19B', 'H193': 'H19C',
                      'H211': 'H21A', 'H212': 'H21B', 'H213': 'H21C',
                      'H221': 'H22A', 'H222': 'H22B',
                      'H231': 'H23A', 'H232': 'H23B',
                      'H241': 'H24A', 'H242': 'H24B',
                      'H261': 'H26A', 'H262': 'H26B', 'H263': 'H26C',
                      'H271': 'H27A', 'H272': 'H27B', 'H273': 'H27C'}
        for n in clol_names :
            atomsel('user 1.0 and resid %d and name %s' % (resid,n)).set('name',clol_names[n])
        A.set('resname','CLOL') # cgenff naming convention
    elif 'BGC' in A.get('resname') :
        print("INFO: Found a beta-glucose!")
        warning = False
        A.set('resname','BGLC')

    # Save in a temporary file
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_extras_', dir=tmp_dir)[1]
    A.set('user', 0.0)
    A.write('pdb',temp)
    name = A.get('resname')[0]

    inp = "whatshappenin"
    if warning:
      print("WARNING! Residue %s isn't explicitly programmed into dabble." % name)
      print("         Please check the atom and residue names in the temp")
      print("         pdb file match those in the topology file before continuing")
      print("         Temp file name: %s" % temp)
      while True :
          if inp in ["Yes","No"]: break
          inp = raw_input("Type Yes to continue, No to skip this residue > ")

    if inp in ["Yes","whatshappenin"] :
        string = '''
   segment %s {
      first none
      last none
      pdb %s
   }
   coordpdb %s %s
        ''' % (name, temp, temp, name)
        file.write(string)
    else :
      os.remove(temp)

    return


# Writes a pdb in order of residues, renumbering the atoms
# accordingly, since psfgen wants each residue sequentially
def write_ordered_pdb(filename, sel, molid=0) :
    f = open(filename, 'w') 
    resids = set( atomsel(sel).get('residue') ) # Much much faster then get resids
    idx = 1
    resnum = 1 # For renumbering residues
    for r in resids :
       # Handle bug where capping groups in same residue as the neighboring amino acid
       # Maestro writes it this way for some reason but it causes problems down the line when
       # psfgen doesn't understand the weird combined residue
       names = set( atomsel('residue %d'% r).get('resname') )
       assert len(names) < 3, "More than 2 residues with same number... currently unhandled. Report a bug"

       if len( names ) > 1 :
          if 'ACE' in names and 'NMA' in names :
              print("ERROR: Both ACE and NMA were given the same residue "
                    "       Check your input structure")
              quit(1)

          if 'ACE' in names :
              names.remove('ACE')
              atomsel('residue %d and resname ACE'% r).set('resid',resnum)
              atomsel('residue %d and resname %s'% (r,names.pop())).set('resid',resnum+1)
              # Handle all of the ACE atoms before the others
              atoms = atomsel('residue %d and resname ACE'% r).get('index')
              atoms.extend( atomsel('residue %d and not resname ACE'% r).get('index') )

          elif 'NMA' in names :
              print("HI ITS A NMA")
              names.remove('NMA')
              atomsel('residue %d and resname %s'% (r,names.pop())).set('resid',resnum)
              atomsel('residue %d and resname NMA'% r).set('resid',resnum+1)
              # Handle all the NMA atoms after the others
              atoms = atomsel('residue %d and not resname NMA'% r).get('index')
              atoms.extend( atomsel('residue %d and resname NMA'% r).get('index') )
          resnum +=2

       else :
           atomsel('residue %d'% r).set('resid',resnum)
           resnum += 1
           atoms = atomsel('residue %d' % r).get('index')

       for i in atoms :
           a = atomsel('index %d' % i)
           
        #  entry="%-6s%5s %5s%-4s%c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n",
           entry='%-6s%5d %-5s%-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n' % ('ATOM',idx,
                 a.get('name')[0], a.get('resname')[0],
                 a.get('chain')[0], a.get('resid')[0],
                 a.get('x')[0], a.get('y')[0], a.get('z')[0],
                 0.0,0.0, a.get('segname')[0], a.get('element')[0] )
           idx += 1
           f.write(entry)
    f.write('END\n')
    atomsel(sel).set('user', 0.0) # Mark as written
    f.close()


# Scans the output psb from psfgen for atoms where the coordinate
# could not be set, because sometimes there is no warning.
def check_psf_output(psf_name) :
    problem_lines = []
    file = open('%s.pdb'% psf_name, 'r')
    for line in iter(file) :
        # Use rsplit since some fields near beginning can mush together
        words = line.rsplit()
        if words[0] == 'ATOM' and words[-4] == '-1.00' :
            problem_lines.append(line[:-1])
    file.close()

    # Print out error messages
    if len(problem_lines) :
        print("\nERROR: Couldn't find the following atoms.\n"
              "       Check if they are present in the original structure.\n"
              "       If they are, check dabble name translation or file a\n"
              "       bug report to Robin.\n")
        for l in problem_lines : 
            print l
    else :
        print("\nINFO: Checked output pdb/psf has all atoms present.\n") 


