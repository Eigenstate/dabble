import molecule
import sys,glob,os,tempfile
from atomsel import *

# Writes a psf file with the given options
def write_psf(psf_name, molid=0, lipid_sel="lipid"):
    """
    Writes a pdb/psf file pair from the current molecule using the CHARMM36
    topology and atom names/types. Interfaces with psfgen by dynamically generating
    the .tcl file that psfgen takes as input. Prompts the user for additional
    topology files and helps with matching atom names that cannot be automatically
    translated to the charmm naming conventions.
    TODO: Someday replace this with ParmEd charmm topology parsing to eliminate
    psfgen

    Args:
      psf_name (str): Prefix for the pdb/psf output files, extension will be appended
      molid (str): the VMD molecule id to write
      lipid_sel (str): atomselect string describing what should count as "lipid"

    Returns:
      topologies (list of str): Topology files that were used in creating the psf
    """

    # Clean up all temp files from previous runs if present
    try :
        os.remove('%s.pdb'% psf_name)
        os.remove('%s.psf'% psf_name)
    except OSError :
        pass

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
    ''' % psf_name
    file.write(string)

    # Specify default parameter sets
    topologies = [ os.environ['DABBLEDIR']+'/charmm_params/top_all36_caps.rtf',
                   os.environ['DABBLEDIR']+'/charmm_params/top_water_ions.rtf',
                   os.environ['DABBLEDIR']+'/charmm_params/top_all36_cgenff.rtf',
                   os.environ['DABBLEDIR']+'/charmm_params/top_all36_prot.rtf',
                   os.environ['DABBLEDIR']+'/charmm_params/top_all36_lipid.rtf',
                   os.environ['DABBLEDIR']+'/charmm_params/top_all36_carb.rtf' ]

    # Ask the user for additional topology files
    sys.stdout.flush()
    print("\nWhich CHARMM topology files should I use?\n"
            "Currently using:")
    for t in topologies :
        print("  - %s" % t.split("/")[-1])

    print("Enter the path to the filename(s) from the current working directory,"
          " separated by a comma, of any other rtf files you wish to use.\n")
    sys.stdout.flush()
    inp = raw_input('> ')
    if inp :
        topologies.extend(inp.split(','))

    file.write('\n')
    for t in topologies :
        file.write('   topology %s\n' % t)

    # Mark all atoms as unsaved with the user field
    atomsel('all',molid=molid).set('user',1.0)

    # Now ions if present, changing the atom names
    if len(atomsel('element Na Cl K',molid=molid)) > 0 :
        temp = write_ion_blocks(file,tmp_dir,molid=molid)
    # End ions

    # Save water 10k molecules at a time
    if len(atomsel('water',molid=molid)) :
        write_water_blocks(file,tmp_dir,molid=molid)
    # End water

    # Now lipid
    if len(atomsel(lipid_sel)) :
        temp=write_lipid_blocks(file,tmp_dir,lipid_sel=lipid_sel,molid=molid)
    # End lipid

    # Here we redefine the protein segments, as if there are breaks in it with
    # capping groups at each end it can still be listed as the same segment by Maestro
    # but since it is not strictly connected should not be the same segment as psfgen 
    # will connect together all atoms in a linear chain.
    # We pick this by resid not residue since sometimes it gets mixed up parsing residues
    # it doesn't know and will split them into multiple residues
    if len(atomsel('protein or resname ACE NMA')) :
        print("Redefining capping groups")
        cap_n = sorted(set(atomsel('resname ACE').get('resid')))
        cap_c = sorted(set(atomsel('resname NMA').get('resid')))
        print cap_n, cap_c

        # Check protein actually has caps and they are done correctly
        if len(cap_c) != len(cap_n) :
          print("\nERROR: There are an uneven number of capping groups on the protein.\n"
                  "       Found %d NMA residues but %d ACE residues\n"
                  "       Please check your input protein preparation.\n" %(len(cap_c),len(cap_n)) )
          quit(1)
        if not len(cap_c) :
          print("\nERROR: There are no capping groups found on the protein.\n"
                  "       Please prepare your input protein with ACE and NMA on the\n"
                  "       terminal residues of each chain segment.\n")
          quit(1);

        segnum=1
        while len(cap_c) > 0 :
            nid = cap_n.pop(0)
            cid = cap_c.pop(0)
            assert nid < cid, "N-terminal residue before C-terminal? NMA resid %d, ACE resid %d"% (nid,cid)
            print("Writing protein chain from resid %d to %d"% (nid,cid))
            sys.stdout.flush()
            segment = atomsel('resid >= %d and resid <= %d and user 1.0'% (nid,cid))
            segment.set('segname','P%s' % segnum)
            segment.set('segid','P%s' % segnum)
            # This is stupid, but to fix the Maestro issue with combined capping group and regualr
            # residues, we have to renumber the residues, save, and reload so that VMD splits the 
            # combined residue into 2
            temp = tempfile.mkstemp(suffix='_P%s.pdb' % segnum, prefix='psf_prot_', dir=tmp_dir)[1]
            segment.write('pdb',temp)
            #write_ordered_pdb(temp, sel='segname P%s'% segnum, molid=molid)
            sys.stdout.flush()

            prot_molid = molecule.load('pdb', temp)
            segment.set('user', 0.0)
            write_protein_blocks(file, tmp_dir=tmp_dir, seg='P%s'%segnum, molid=prot_molid, topologies=topologies)
            molecule.delete(prot_molid)
            #os.remove(temp)
            segnum += 1
    else :
      print("\nINFO: Didn't find any protein.\n")
    # End protein

    # Check if there is anything else and let the user know about it
    leftovers = atomsel('user 1.0', molid=molid)
    if len(leftovers) > 0 :
        print("Found extra ligands: %s" % set(leftovers.get('resname')))
    for l in set( leftovers.get('resid') ) :
        write_ligand_blocks(file, tmp_dir, resid=l, topologies=topologies, molid=molid)

    # Write the output files and run
    string='''
    writepsf x-plor cmap ${output}.psf
    writepdb ${output}.pdb'''
    file.write(string)
    file.close()

    from VMD import evaltcl
    evaltcl('play %s' % filename)
    check_psf_output(psf_name)

    return topologies

# Writes the temporary protein PDB file with correct atom names for psfgen
# The molid is of a molecule containing only the protein chain to be written
def write_protein_blocks(file,tmp_dir,seg,molid,topologies):
    # Put molid on top to simplify atom selections
    old_top=molecule.get_top()
    molecule.set_top(molid)
    print("Writing protein file\n")

    # Renumber residues starting from 1
    residues = set( atomsel('all').get('residue') )
    resnum = 1
    while len(residues) :
        atomsel('residue %s'% residues.pop()).set('resid',resnum)
        resnum += 1
    #T.set('resid',[r+1 for r in T.get('residue')])

    # Check the protein has hydrogens
    if len(atomsel('element H')) == 0:
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
        atomsel('resname ACE and name %s'% n).set('name',ace_names[n])

    # Terminal residue NMA
    nma_names = {'HN'  :'HNT', 'H'   :'HNT', 
                 'CT3' :'CAT', 'CA'  :'CAT', 'CH3' :'CAT',
                 'HA1' :'HT1', 'HA2' :'HT2', 'HA3' :'HT3',
                 '1HA' :'HT1', '2HA' :'HT2', '3HA' :'HT3', 
                 'HH31':'HT1', 'HH32':'HT2', 'HH33':'HT3' }
    for n in nma_names :
        atomsel('resname NMA and name %s' % n).set('name',nma_names[n])

    patches = ''

    # Disulfide briges
    indices = set( atomsel('name SG and resname CYX').get('index') )
    indices.update( atomsel('name SG and resname CYS and not same residue as name HG').get('index') )
    while len(indices) > 0 :
        idx1 = indices.pop()
        matches = set( atomsel('name SG and (resname CYX or (resname CYS and not same residue as name HG)) and not index %d and within 2.5 of index %d' % (idx1, idx1) ))

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
    atomsel('resname HID').set('resname','HSD')
    atomsel('resname HIE').set('resname','HSE')
    atomsel('resname HIP').set('resname','HSP') 

    # Determine protonation states of those just called HIS based on atom names
    atomsel('resname HIS and same residue as name HE2').set('resname','HSE')
    atomsel('resname HIS and same residue as name HD1').set('resname','HSD')
    # NOTE: This atomsel MUST come after the previous two!
    atomsel('resname HIS and same residue as name HE2 and same residue as name HD1').set('resname','HSP')
    
    # Isoleucine
    iso_names = { 'CD1' :  'CD',
                  'HD11': 'HD1', 'HD12':'HD2', 'HD13':'HD3',
                  'HG12':'HG11', 'HG13':'HG12' }
    for n in iso_names :
        atomsel('resname ILE and name %s' %n).set('name',iso_names[n])
                  
    # Glutamine check all residues for protonation
    # Can't use dictionary for names since two of them must be swapped
    # Must loop by resid, not residue, to get correct PATCH statement index
    t = set( atomsel('resname GLU GLH GLUP').get('resid') )
    for resid in t:
        atomsel('resid %s' % resid).set('resname','GLU')
        if "HE1" in atomsel('resid %s' % resid).get('name') :
            atomsel('resid %s and name HE1'% resid).set('name','HE2')
            atomsel('resid %s and name OE1'% resid).set('name','temp')
            atomsel('resid %s and name OE2'% resid).set('name','OE1')
            atomsel('resid %s and name temp'% resid).set('name','OE2')
            patches += 'patch GLUP %s:%d\n' % (seg, resid)
        elif "HE2" in atomsel('resid %s' % resid).get('name') :
            patches += 'patch GLUP %s:%d\n' % (seg, resid)
        elif "HXT" in atomsel('resid %s' % resid).get('name') :
            atomsel('resid %s and name HXT' % resid).set('name','HE2')
            patches += 'patch GLUP %s:%d\n' % (seg, resid)

    # Aspartate check each residue to see if it's protonated
    # Can't use dictionary for names since two of them must be swapped
    # Must loop by resid, not residue, to get correct PATCH statement index
    t = set( atomsel('resname ASP ASH ASPP').get('resid') )
    for resid in t :
        atomsel('resid %s'% resid).set('resname','ASP')
        if "HD1" in atomsel('resid %s' % resid).get('name') :
            atomsel('resid %s and name HD1' % resid).set('name','HD2')
            atomsel('resid %s and name OD1' % resid).set('name','temp')
            atomsel('resid %s and name OD2' % resid).set('name','OD1')
            atomsel('resid %s and name temp' % resid).set('name','OD2')
            patches += 'patch ASPP %s:%d\n' % (seg,resid)
        if "HD2" in atomsel('resid %s' % resid).get('name') :
            patches += 'patch ASPP %s:%d\n' % (seg,resid)

    # Methionine hydrogen names
    atomsel('name H2 H1 and resname MET').set('name','HN')

    # Serine and cysteine hydrogen names
    atomsel('name HG and resname SER CYS').set('name','HG1')

    # Hydrogens can be somewhat residue-dependent, but only check hydrogens on known amino acid
    # residues so that the names on nonstandard amino acids are never changed
    acids = 'ACE ALA ARG ASN ASP CYS GLN GLU GLY HIS HSP HSE HSD ILE LEU LYS MET NMA PHE PRO SER THR TRP TYR VAL'
    h_names = {'H'  : 'HN', 'H2' :'HN',
               'HA2':'HA1',
               'HA3':'HA2', 'HG2':'HG1',
               'HG3':'HG2'}
    for h in h_names :
        atomsel('resname %s and name %s' % (acids,h)).set('name',h_names[h])

    # These two statements must execute in this specific order
    atomsel('resname %s and name HB2 and not resname ALA'% acids).set('name','HB1')
    atomsel('resname %s and name HB3 and not resname ALA'% acids).set('name','HB2')

    # Some residue-specific stuff
    atomsel('resname %s and name HD2 and not (resname TYR ILE HSP HSE HSD ASP PHE)'
            % acids).set('name','HD1')
    atomsel('resname %s and name HD3 and not (resname TYR ILE HIS HSE HSD ASP PHE)' 
            % acids).set('name','HD2')
    atomsel('resname %s and name HE2 and not (resname TYR MET TRP HSP HSE HSD GLU PHE)'
            % acids).set('name','HE1')
    atomsel('resname %s and name HE3 and not (resname TYR MET TRP HSP HSE HSD PHE)' 
            % acids).set('name','HE2')

    # Final chance to fix unrecognized atom names for non-protein residues
    others = set( atomsel('not resname %s'% acids).get('resname') )
    while len(others) :
        res = others.pop()
        # Check if this residue is not found, and offer to rename it
        if not find_residue_in_rtf(topologies=topologies,resname=res,segname=seg,molid=molid) :
            print("\nERROR: Residue name %s wasn't found in any input topology.\n"
                    "       Would you like to rename it?\n" % res)
            newname = raw_input("New residue name or CTRL+D to quit > ")
            atomsel('resname %s'% res).set('resname','%s'% newname)
            others.add(newname) # Check that new name is present

    # Now protein
    filename = tmp_dir + '/psf_protein_%s.pdb'% seg
    write_ordered_pdb(filename, sel='all', molid=molid)
    print("Wrote %d atoms to the protein segment %s" % (len(atomsel('all' )),seg))
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

    # Detect what lipid kind this is. HACK for now, it really should figure it out,
    # but this is what works
    all = atomsel('(%s) and user 1.0' % lipid_sel, molid=molid)
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
    # These selections must be done in order because there is a swap
    atomsel('(%s) and name C11' % lipid_sel).set('name','t12')
    atomsel('(%s) and name C15' % lipid_sel).set('name','C11')
    atomsel('(%s) and name C14' % lipid_sel).set('name','C15')
    atomsel('(%s) and name C12' % lipid_sel).set('name','C14')
    atomsel('(%s) and name t12' % lipid_sel).set('name','C12')

    popc_names = {'H31':'H13A', 'H32':'H13B', 'H33':'H13C',
                  'H41':'H15A', 'H42':'H15B', 'H43':'H15C',
                  'H21':'H14A', 'H22':'H14B', 'H23':'H14C',
                  'H51':'H11A', 'H52':'H11B',
                  'H11':'H12A', 'H12':'H12B',
                 # Phosphate and its oxygens
                  'P1' :'P'   , 'O1' :'O12' , 'O2' :'O11',
                  'O3' :'O13' , 'O4' :'O14' }
    # Fortunately all the other atom names are correct.
    for n in popc_names :
        atomsel('(%s) and name %s' % (lipid_sel,n)).set('name',popc_names[n])

    # Write temporary lipid pdb
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_lipid_', dir=tmp_dir)[1]
    all.set('user',0.0)
    all.write('pdb',temp)

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


def write_ligand_blocks(file, tmp_dir, resid, topologies, molid=0) :
    A = atomsel('user 1.0 and resid %d' % resid)

    # Adjust residue name if known ligand
    warning = True
    if 'CLR' in A.get('resname') :
        print("INFO: Found a cholesterol!")
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
    #elif 'BGC' in A.get('resname') :
    #    print("INFO: Found a beta-glucose!")
    #    A.set('resname','BGLC')
    else : # Check it's recognized, if not give the user a chance to rename
         res = A.get('resname')[0]
         seg = A.get('segname')[0]
         while not find_residue_in_rtf(topologies=topologies,resname=res,segname=seg,molid=molid) :
             print("\nERROR: Residue name %s wasn't found in any input topology.\n"
                     "       Would you like to rename it?\n" % res)
             newname = raw_input("New residue name or CTRL+D to quit > ")
             atomsel('resid %s and resname %s'% (seg,res)).set('resname','%s'% newname)

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
# could not be set, because sometimes there is no warning in the output text.
def check_psf_output(psf_name) :
    # Check file was written at all
    if not os.path.isfile('%s.pdb'% psf_name) :
        print("\nERROR: psf file failed to write.\n"
                "       Please see log above.\n")
        quit(1)

    problem_lines = []
    file = open('%s.pdb'% psf_name, 'r')
    for line in iter(file) :
        # Use rsplit since some fields near beginning can mush together
        # 3rd field from end not fourth since element name will be absent too
        words = line.rsplit()
        if words[0] == 'ATOM' and words[-3] == '-1.00' :
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
        quit(1)
    else :
        print("\nINFO: Checked output pdb/psf has all atoms present and correct.\n") 


# Scans the input text for the residue with the given name.
# Once found, pulls out all atom names that comprise that residue and returns them as a set.
# If the residue isn't present in the input text, returns an empty set.
def get_atoms_from_rtf(text,resname) :
    atoms = []
    found = False
    for i in range(len(text)) :
        words = text[i].split()
        if not len(words)             : continue
        if not found and words[0]=='RESI' \
           and words[1]==resname      : found = True
        elif found and words[0]=='ATOM' : 
            atoms.append(words[1])
        elif found and words[0]=='RESI' : break
    return set(atoms)

# Scan the topology files to find this residue name
# Pulls out the atoms involved and checks that they are all present
# in the input coordinates. Allows name correction if not.
# Returns False if the residue name cannot be found
def find_residue_in_rtf(topologies,resname,segname,molid) :
    names = ()
    found = False
    for t in topologies :
        topfile = open(t, 'r')
        topo_atoms = get_atoms_from_rtf(text=topfile.readlines(),resname=resname)
        if len(topo_atoms) : break # Use first definition found of this residue
        topfile.close()
    if not len(topo_atoms) : return False
    print("INFO: Successfully found residue %s in input topologies"% resname)

    # Match up atoms with python sets
    pdb_atoms = set(atomsel('segname %s and resname %s'% (segname,resname)).get('name'))
    pdb_only  = pdb_atoms - topo_atoms
    topo_only = topo_atoms - pdb_atoms

    # If uneven number of atoms, there are missing or additional atoms
    if len(pdb_atoms) > len(topo_atoms) :
        print("\nERROR: Cannot process modified residue %s.\n"
                "       There are %d extra atoms in the input structure that are \n"
                "       undefined in the topology file. The following atoms could not \n"
                "       be matched and may either be misnamed, or additional atoms.\n"
                "       Please check your input." % (resname,len(pdb_atoms)-len(topo_atoms)))
        print(  "       [ %s ]\n" % ' '.join(map(str,pdb_only)) )
        print(  "       Cannot continue.\n")
        quit(1);
    if len(topo_atoms) > len(pdb_atoms) :
        print("\nERROR: Cannot process modified residue %s.\n"
                "       There are %d missing atoms in the input structure that are \n"
                "       defined in the topology file. The following atoms could not \n"
                "       be matched and may either be misnamed or deleted atoms.\n"
                "       Please check your input." % (resname,len(topo_atoms)-len(pdb_atoms)))
        print(  "       [ %s ]\n" % ' '.join(map(str,topo_only)) )
        print(  "       Cannot continue.\n")
        quit(1);

    # Offer to rename atoms that couldn't be matched to the topology
    if len(pdb_only) :
        print("\nWARNING: Having some trouble with modified residue %s.\n"
                "         The following atom names cannot be matched up to the input\n"
                "         topologies. They are probably misnamed.\n" % resname)
        print(  "         To help you, here are the atom names that should be present\n"
                "         according to the topology but were not found:\n")
        print(  "         [ %s ]\n" % ' '.join(map(str,topo_only)) )
        print(  "         Please enter a valid name for each atom as it appears or CTRL+D to quit..\n")
        for r in pdb_only :
            print("Unmatched topology names: [ %s ]" % ' '.join(map(str,topo_only)) )
            newname = raw_input("  %s  -> "% r)
            while newname not in topo_only:
                print("'%s' is not an available name in the topology. Please try again.\n" % newname)
                newname = raw_input("  %s  -> "% r)
            atomsel('segname %s and resname %s and name %s'% (segname,resname,r)).set('name',newname)
            pdb_atoms = set(atomsel('segname %s and resname %s'% (segname,resname)).get('name'))
            topo_only = topo_atoms-pdb_atoms

        # Recurse to check that everything is assigned correctly
        find_residue_in_rtf(topologies,resname,segname,molid)
    print("INFO: Matched up all atom names for residue %s\n"% resname)
    return True


def write_amber(psf_name, prmtop_name, topologies, final_molid) :
    """
    Runs the chamber command of ParmEd to produce AMBER format input files.

    Args:
      psf_name (str): The file name of the psf and pdb files to convert. Should
          just include the prefix, not the .psf or .pdb extension.
      prmtop_name (str): The prefix for the output .prmtop and .inpcrd files.
          The extension will be appended.
      topologies (list of str): CHARMM format topology files that define the system,
          including default list and any user-defined ones.
      final_molid (int): VMD molid of loaded molecule containing the final system.
          Used for box information as it is not written by psfgen to the pdb/psf.

    Returns:
      True if successful
    """ 

    # Ask the user for additional parameter files
    parameters = [ os.environ['DABBLEDIR']+'/charmm_params/toppar_water_ions.str',
                   os.environ['DABBLEDIR']+'/charmm_params/par_all36_cgenff.prm',
                   os.environ['DABBLEDIR']+'/charmm_params/par_all36_prot.prm',
                   os.environ['DABBLEDIR']+'/charmm_params/par_all36_lipid.prm',
                   os.environ['DABBLEDIR']+'/charmm_params/par_all36_carb.prm' ]

    sys.stdout.flush()
    print("\nWhich CHARMM parameter files should I use?\n"
            "Currently using:")
    for p in parameters:
        print("  - %s" % p.split("/")[-1])
    print("Enter the path to the filename(s) from the current working directory,"
          " separated by a comma, of any other rtf files you wish to use.\n")
    sys.stdout.flush()
    inp = raw_input('> ')
    if inp :
        parameters.extend(inp.split(','))

    # Begin assembling chamber input string
    args = "-crd %s.pdb -psf %s.psf" % (psf_name,psf_name)

    # Add topology and parameter arguments
    for t in topologies :
      args += ' -top %s' % t
    for p in parameters :
      if 'toppar' in p : args += ' -toppar %s' %p
      else :             args += ' -param %s' % p


    # Add box information since it is not in the pdb
    box = molecule.get_periodic(molid=final_molid)
    args += " -box %f,%f,%f" % (box['a'],box['b'],box['c']) 

    from ParmedTools import chamber, parmout
    from chemistry.amber import AmberParm
    print("\nINFO: Running chamber. This may take a while...")
    sys.stdout.flush()
    parm = AmberParm()
    action = chamber(parm,args)
    print action
    action.execute()
    print("\nINFO: Ran chamber")
    write = parmout(action.parm, "%s.prmtop %s.inpcrd"%(prmtop_name,prmtop_name))
    write.execute()
    print("\nINFO: Wrote output prmtop and inpcrd")
    return True



