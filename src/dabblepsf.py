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
       - top_all36_lipid.rtf """)
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

    # Save the protein with correct atom names
    chains = set( atomsel('user 1.0 and (protein or (resname ACE NMA))',molid=molid).get('chain') )
    for ch in chains : 
        print("Writing protein chain %s" % ch)
        write_protein_blocks(file,chain=ch,tmp_dir=tmp_dir,molid=molid)
    # End protein

    # Check if there is anything else and let the user know about it
    print("Found extra ligands: %s" % set( atomsel('user 1.0',molid=molid).get('resname') ))
    leftovers = set( atomsel('user 1.0',molid=molid).get('residue') )
    for l in leftovers :
        write_ligand_blocks(file, tmp_dir, residue=l, molid=molid)

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
def write_protein_blocks(file,tmp_dir,chain,molid=0):
    # Put molid on top to simplify atom selections
    old_top=molecule.get_top()
    molecule.set_top(molid)

    # Renumber residues starting from 1
    T = atomsel('(protein or resname ACE NMA) and chain %s' % chain)
    T.set('resid',[r+1 for r in T.get('residue')])

    # Check the protein has hydrogens
    if len(atomsel('protein and element H')) == 0:
        print("\n\nERROR: There are no hydrogens on your protein")
        print("       You need to add them before parameterizing")
        print("       because psfgen cannot be trusted to do it correctly.\n")
        quit(1)

    # Terminal residue ACE
    atomsel('chain %s and resname ACE and name CH3' % chain).set('name','CAY')
    atomsel('chain %s and resname ACE and name HH31'% chain).set('name','HY1') # H could be named HH31 or H1
    atomsel('chain %s and resname ACE and name HH32'% chain).set('name','HY2')
    atomsel('chain %s and resname ACE and name HH33'% chain).set('name','HY3')
    atomsel('chain %s and resname ACE and name H1'% chain).set('name','HY1')
    atomsel('chain %s and resname ACE and name H2'% chain).set('name','HY2')
    atomsel('chain %s and resname ACE and name H3'% chain).set('name','HY3')
    atomsel('chain %s and resname ACE and name 1H'% chain).set('name','HY1')
    atomsel('chain %s and resname ACE and name 2H'% chain).set('name','HY2')
    atomsel('chain %s and resname ACE and name 3H'% chain).set('name','HY3')
    atomsel('chain %s and resname ACE and name C'% chain).set('name','CY')
    atomsel('chain %s and resname ACE and name O'% chain).set('name','OY')

    # Terminal residue NMA
    atomsel('chain %s and resname NMA and name HN'% chain).set('name','HNT')
    atomsel('chain %s and resname NMA and name H'% chain).set('name','HNT')
    atomsel('chain %s and resname NMA and name N'% chain).set('name','NT')
    atomsel('chain %s and resname NMA and name CA'% chain).set('name','CAT')
    atomsel('chain %s and resname NMA and name HA1'% chain).set('name','HT1')
    atomsel('chain %s and resname NMA and name HA2'% chain).set('name','HT2')
    atomsel('chain %s and resname NMA and name HA3'% chain).set('name','HT3')

    patches = ''
    # NOTE: No patch for ACE or NMA  since this is done in the termini rtf file as a regular residue

    # Disulfide briges
    indices = set( atomsel('chain %s and name SG and resname CYX'% chain).get('index') )
    while len(indices) > 0 :
        idx1 = indices.pop()
        matches = set( atomsel('chain %s and name SG and resname CYX and not index %d and within 2.5 of index %d' % (chain, idx1, idx1) ))

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
        patches += 'patch DISU P:%d P:%d\n' % (res1,res2)

    # Rename CYX -> CYS for naming conventions
    atomsel('resname CYX').set('resname','CYS')

    # Histidine naming convention
    atomsel('chain %s and resname HID'% chain).set('resname','HSD')
    atomsel('chain %s and resname HIE'% chain).set('resname','HSE')
    atomsel('chain %s and resname HIP'% chain).set('resname','HSP') 

    # Determine protonation states of those just called HIS based on atom names
    atomsel('chain %s and resname HIS and same residue as name HE2'% chain).set('resname','HSE')
    atomsel('chain %s and resname HIS and same residue as name HD1'% chain).set('resname','HSD') # NOTE: This atomsel MUST come after the previous two!
    atomsel('chain %s and resname HIS and same residue as name HE2 and same residue as name HD1' % chain).set('resname','HSP')
    
    # Isoleucine
    atomsel('chain %s and resname ILE and name CD1'% chain).set('name','CD')
    atomsel('chain %s and resname ILE and name HD11'% chain).set('name','HD1')
    atomsel('chain %s and resname ILE and name HD12'% chain).set('name','HD2')
    atomsel('chain %s and resname ILE and name HD13'% chain).set('name','HD3')
    atomsel('chain %s and resname ILE and name HG12'% chain).set('name','HG11')
    atomsel('chain %s and resname ILE and name HG13'% chain).set('name','HG12')

    # Glutamine check all residues for protonated one
    t=atomsel('chain %s and (resname GLU or resname GLH)'% chain).get('residue')
    for residue in t:
        if "HE1" in atomsel('residue %s' % residue).get('name') :
            atomsel('chain %s and residue %s and name HE1' % (chain,residue)).set('name','HE2')
            atomsel('chain %s and residue %s and name OE1' % (chain,residue)).set('name','temp')
            atomsel('chain %s and residue %s and name OE2' % (chain,residue)).set('name','OE1')
            atomsel('chain %s and residue %s and name temp' % (chain,residue)).set('name','OE2')
            atomsel('chain %s and residue %s').set('resname','GLUP')
        elif "HE2" in atomsel('chain %s and residue %s' % (chain,residue)).get('name') :
            atomsel('chain %s and residue %s' % (chain,residue)).set('resname','GLUP')
        # TODO: Is this correct? Atom name messups
        elif "HXT" in atomsel('chain %s and residue %s' % (chain,residue)).get('name') :
            atomsel('chain %s and residue %s' % (chain,residue)).set('resname','GLUP')
            atomsel('chain %s and residue %s and name HXT' % (chain,residue)).set('name','HE2')
        else :
            atomsel('chain %s and residue %s' % (chain,residue)).set('resname','GLU')

    # Aspartate check each residue to see if it's protonated
    t=atomsel('resname ASP or resname ASH').get('residue')
    for residue in t :
        if "HD1" in atomsel('chain %s and residue %s' % (chain,residue)).get('name') :
            atomsel('chain %s and residue %s and name HD1' % (chain,residue)).set('name','HD2')
            atomsel('chain %s and residue %s and name OD1' % (chain,residue)).set('name','temp')
            atomsel('chain %s and residue %s and name OD2' % (chain,residue)).set('name','OD1')
            atomsel('chain %s and residue %s and name temp' % (chain,residue)).set('name','OD2')
            atomsel('chain %s and residue %s').set('resname','ASPP')
        if "HD2" in atomsel('chain %s and residue %s' % (chain,residue)).get('name') :
            atomsel('chain %s and residue %s' % (chain,residue)).set('resname','ASPP')
        else :
            atomsel('chain %s and residue %s' % (chain,residue)).set('resname','ASP')

    for resid in set(atomsel('chain %s and resname ASPP' % chain).get('resid')) :
        patches += 'patch ASPP P:%d\n' % resid
    atomsel('resname ASPP').set('resname','ASP')

    # Hydrogens
    atomsel('chain %s and (resname ACE or protein) and name H' % chain).set('name','HN')
    atomsel('chain %s and (resname ACE or protein) and name HA2' % chain).set('name','HA1')
    atomsel('chain %s and (resname ACE or protein) and name HA3' % chain).set('name','HA2')
    atomsel('chain %s and (resname ACE or protein) and name HG2' % chain).set('name','HG1')
    atomsel('chain %s and (resname ACE or protein) and name HG3' % chain).set('name','HG2')
    atomsel('chain %s and (resname ACE or protein) and name HB2 and not resname ALA' % chain).set('name','HB1')
    atomsel('chain %s and (resname ACE or protein) and name HB3 and not resname ALA' % chain).set('name','HB2')
    atomsel('chain %s and (resname ACE or protein) and name HD2 and not (resname ILE HSP HSE HSD ASP PHE)' % chain).set('name','HD1')
    atomsel('chain %s and (resname ACE or protein) and name HD3 and not (resname ILE HIS HSE HSD ASP PHE)' % chain).set('name','HD2')
    atomsel('chain %s and (resname ACE or protein) and name HE2 and not (resname TRP HSP HSE HSD GLUP PHE)' % chain).set('name','HE1')
    atomsel('chain %s and (resname ACE or protein) and name HE3 and not (resname TRP HSP HSE HSD PHE)' % chain).set('name','HE2')
    atomsel('chain %s and (resname ACE or protein) and name HG and resname SER CYS' % chain).set('name','HG1')

    # Now protein
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_protein_', dir=tmp_dir)[1]
    print("Writing temporary protein file...")
    write_ordered_pdb(temp, sel='chain %s and (resname ACE NMA or protein)' % chain, molid=molid)
    print("Wrote %d atoms to the protein chain %s" % (len(atomsel('chain %s and resname ACE or protein' % chain)),chain))
    molecule.set_top(old_top)

    # Now write to psfgen input file
    string='''
      set protnam %s
      segment P%s {
        pdb $protnam
      } 
    ''' % (temp,chain)
    file.write(string)
    file.write(patches)
    file.write('coordpdb $protnam P%s' % chain) 

    return temp

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


def write_ligand_blocks(file, tmp_dir, residue, molid=0) :
    A = atomsel('residue %d' % residue)

    # Adjust residue name if known ligand
    if 'CLR' in A.get('resname') :
        A.set('resname','CLOL') # cgenff naming convention

    # Save in a temporary file
    temp = tempfile.mkstemp(suffix='.pdb', prefix='psf_extras_', dir=tmp_dir)[1]
    A.set('user', 0.0)
    A.write('pdb',temp)
    name = A.get('resname')[0]

    string = '''
      segment %s {
        pdb %s
      }
      coordpdb %s %s
    ''' % (name, temp, temp, name)
    file.write(string)
    return temp


# Writes a pdb in order of residues, renumbering the atoms
# accordingly, since psfgen wants each residue sequentially
def write_ordered_pdb(filename, sel, molid=0) :
    f = open(filename, 'w')
    atomsel(sel).set('user', 0.0) # Mark as written

    resids = set( atomsel(sel).get('residue') ) # Much much faster then get resids
    idx = 1
    for r in resids :
       atoms = atomsel('residue %d' % r).get('index')
       for i in atoms :
           a = atomsel('index %d' % i)
           
        #  entry="%-6s%5s %5s%-4s%c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n",
           entry='%-6s%5d %-4s%c%-4s%c%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n' % ('ATOM',idx,
                 a.get('name')[0],' ', a.get('resname')[0],
                 a.get('chain')[0], a.get('resid')[0],' ',
                 a.get('x')[0], a.get('y')[0], a.get('z')[0],
                 0.0,0.0, a.get('segname')[0], a.get('element')[0] )
           idx += 1
           f.write(entry)
    f.close()

