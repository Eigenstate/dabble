# dabblelib.py
# By Robin M. Betz
# 21 January 2015
# Loosely based off saute3.py from DESRES

import os
import tempfile

from numpy import *

import vmd
from atomsel import atomsel

# Yes, these are different modules. With different functions. Only 2nd one currently needed
# VMD is the worst.
#from Molecule import * # http://www.ks.uiuc.edu/Research/vmd/current/ug/node188.html
#from molecule import * # http://www.ks.uiuc.edu/Research/vmd/current/ug/node182.html
import molecule


############################################################
#                       CONSTANTS                          #
############################################################


__MEMBRANE_HYDROPHOBIC_THICKNESS = 30.0
__MEMBRANE_FULL_THICKNESS = 50.0

__1M_SALT_IONS_PER_WATER = 0.018

random.seed(2015)

############################################################
#             STRUCTURE ANALYSIS FUNCTIONS                 #
############################################################

def orient_solute(molid, z_move, z_rotation, opm_pdb=None, opm_align='protein and backbone'):

    # Check that OPM and alignment aren't both specified
    if opm_pdb is not None and (z_move is not 0 or z_rotation is not 0) :
        print("ERROR: Cannot specify an OPM pdb and manual orientation information");
        quit(1);

    if opm_pdb is not None :
        opm = molecule.load('pdb',opm_pdb)
        T=atomsel('protein and backbone', molid=molid).fit(atomsel(opm_align,molid=opm))
        atomsel('all', molid=molid).move(T)
        molecule.delete(opm)
        return molid

    import trans, math
    if z_move is not None :
        print("Moving z!")
        atomsel('all',molid=molid).moveby((0,0,z_move))
        if z_rotation is None:
          return molid

    if z_rotation is not None:
        print("Rotating!")
        trans.resetview(molid) # View affect rotation matrix, now it's I
        # This is negative because we want membrane flat along the z-axis,
        # and OPM lists the membrane rotation relative to the protein
        theta = math.radians(-1*z_rotation)
        # Rotation matrix in row order with 4th dimension just from I
        rotmat = [ math.cos(theta), -1*math.sin(theta), 0, 0,
                   math.sin(theta),   math.cos(theta),  0, 0,
                         0        ,         0,          1, 0,
                         0        ,         0,          0, 1 ]
        trans.set_rotation(molid, rotmat)
        return molid

    molid=center_system(molid=molid, center_z=False)

    return molid


def get_net_charge(sel,molid):
    """
    Gets the net charge of an atom selection, using the charge
    field of the data.

    Args:
      sel (str): VMD atom selection to compute the charge of
      molid (int): VMD molecule id to select within

    Returns:
      (int): The rounded net charge of the selection

    Throws:
      AssertionException: If all charges are zero
    """

    charge = array(atomsel(sel).get('charge'))
    if charge.size == 0:
        return 0
    if all(charge == 0) :
        print("\nWARNING: All charges in selection are zero. Check the input file has formal charges defined!\nSelection was:\n%s\n"%sel)
        print set(charge)
    #assert all (charge == 0.), 'All charges in selection "%s" are zero! Check your input file has charges'%sel
    net_charge = sum(charge)
    rslt = round(net_charge)
    assert abs(rslt - net_charge) < .01, 'Total charge is not integral (within a tolerance of .01). Check your input file.'
    return int(rslt)

def get_system_net_charge(molid):
    """
    Gets the net charge of the entire system.
    What is in the system is defined by the beta field, as atoms that won't
    be written have beta 0.

    Args:
      molid (int): VMD molecule id to compute the charge of

    Returns:
      (int): The net charge of the molecule
    """

    return get_net_charge(sel='beta 1', molid=molid)


def diameter(v, chunkmem=30e6):
    import sys
    v = require(v, dtype=float32)
    if v.size == 0:
        return 0
    nper = int(chunkmem / (4*v.size))
    D = array([inner(x,x) for x in v])
    i = 0
    d = 0
    while i < len(v):
#        sys.stderr.write('%d..%d\n' % (i, i+nper))
        M = inner(-2*v, v[i:i+nper])
        M += D[i:i+nper]
        M += D[:,None]
        nd = M.max()
        if nd > d:
#            sys.stderr.write('%f -> %f\n' % (d, nd))
            d = nd
        i += nper
    return sqrt(d)


def solute_xy_diameter(solute_sel):
    a = atomsel(solute_sel)
    return diameter(transpose([a.get('x'), a.get('y')]))

# Accepts either a molecule id or file to load. 
def get_solute_sel(molid=0) :
    chains = set(atomsel('all', molid=molid).get('chain'))
    ch = chains.pop() # have to handle first separately because of or
    sel = '(chain %s and resid ' % ch + \
          ' '.join(map(str,set(atomsel('chain %s'% ch, molid=molid).get('resid')))) + ')'
    for ch in chains :
        sel += 'or (chain %s and resid ' % ch + \
               ' '.join(map(str,set(atomsel('chain %s'% ch, molid=molid).get('resid')))) + ')'
    return sel


def get_cell_size(mem_buf,
                  wat_buf,
                  solute_sel,
                  molid=None,
                  filename=None,
                  zh_mem_full=__MEMBRANE_FULL_THICKNESS / 2.0,
                  zh_mem_hyd=__MEMBRANE_HYDROPHOBIC_THICKNESS / 2.0):

    if filename is not None :
        top = molecule.get_top()    
        tmp_top = molecule.read(-1, 'mae', filename)

    if molid is None:
        molid = molecule.get_top()

    solute_z = atomsel(solute_sel, molid=molid).get('z') + [-zh_mem_full, zh_mem_full] # add in dummy values for the membrane boundaries in case the protein is peripheral
    z_min, z_max = min(solute_z), max(solute_z)
    dxy_sol = solute_xy_diameter(solute_sel)
    dxy_tm  = solute_xy_diameter('(%s) and z > %f and z < %f' % (solute_sel, -zh_mem_hyd, zh_mem_hyd))
    xy_size = max(dxy_tm + mem_buf, dxy_sol + wat_buf)
    dz_full = z_max - z_min
    z_size =  dz_full + wat_buf

# Cleanup temporary file, if read in
    if filename is not None :
        molecule.delete(tmp_top)
        if top != -1: molecule.set_top(top)

    return xy_size, z_size, dxy_sol, dxy_tm, dz_full


def get_num_salt_ions_needed(conc,
                             water_sel = 'water and element O',
                             cation='Na',
                             anion='Cl'):
    cations = atomsel_remaining('element %s' % cation)
    anions = atomsel_remaining('element %s' % anion)
    num_cations = len(cations)
    num_anions = len(anions)
    molid = molecule.get_top()
    try:
        abs(get_net_charge(str(cations),molid)-num_cations) > 0.01
    except:
        # Check for bonded cations
        nonbonded_cation_index=[atomsel('index %d' %x).get('index')[0] for x in nonzero(array(map(len,atomsel('element %s' %cation).bonds)) == 0)[0]]
        if len(nonbonded_cation_index)==0:
            cations=atomsel('none')
        else:
            cations=atomsel_remaining('index %s' %' '.join(map(str,nonbonded_cation_index)))
        num_cations = len(cations)
        assert abs(get_net_charge(str(cations),molid)-num_cations) < 0.01, 'num cations and net cation charge are not equal'
    try:
        abs(get_net_charge(str(anions),molid)+num_anions) > 0.01
    except:
        # Check for bonded anions
        non_bonded_anion_index=[atomsel('index %d' %x).get('index')[0] for x in nonzero(array(map(len,atomsel_remaining('element %s' %anion).bonds)) == 0)[0]]
        if len(non_bonded_anion_index)==0:
            anions=atomsel('none')
        else:
            anions=atomsel_remaining('index %s' %' '.join(map(str,nonbonded_anion_index)))
        num_anions = len(anions)
        assert abs(get_net_charge(str(anions),molid)+num_anions) < 0.01, 'num anions and abs anion charge are not equal'
    num_waters = num_atoms_remaining(water_sel)
    num_for_conc = int(round(__1M_SALT_IONS_PER_WATER * num_waters * conc))
    pos_ions_needed = num_for_conc - num_cations
    neg_ions_needed = num_for_conc - num_anions
    system_charge = get_system_net_charge(molid)

    new_system_charge = system_charge + num_anions - num_cations
    to_neutralize = abs(new_system_charge)
    if new_system_charge > 0:
        if to_neutralize > pos_ions_needed:
            neg_ions_needed += to_neutralize - pos_ions_needed
            pos_ions_needed = 0
        else:
            pos_ions_needed -= to_neutralize
    else:
        if to_neutralize > neg_ions_needed:
            pos_ions_needed += to_neutralize - neg_ions_needed
            neg_ions_needed = 0
        neg_ions_needed -= to_neutralize
        
    total_cations = num_cations + pos_ions_needed
    total_anions = num_anions + neg_ions_needed
    
    cation_conc = (float(total_cations) / num_waters) / __1M_SALT_IONS_PER_WATER # volume estimate from prev waters
    anion_conc = (float(total_anions) / num_waters) / __1M_SALT_IONS_PER_WATER
    num_waters -= pos_ions_needed + neg_ions_needed

    return (pos_ions_needed,
            neg_ions_needed,
            num_waters,
            total_cations,
            total_anions,
            cation_conc,
            anion_conc)


def lipid_composition(lipid_sel):
    def leaflet(leaflet_sel):
        sel = atomsel_remaining('not element H C and (%s) and (%s)' % (lipid_sel, leaflet_sel))
        resnames = set(sel.get('resname'))
        d = dict([(s, len(set(atomsel_remaining('not element H C and resname %s and (%s) and (%s)' % (s, lipid_sel, leaflet_sel)).get('fragment')))) for s in resnames])
        return d
    inner, outer = leaflet('z < 0'), leaflet('not (z < 0)')
    return inner, outer


############################################################
#              MAE FILE MANIPULATION FUNCTIONS             #
############################################################

# Loads the input file, either pdb or mae, returns molecule id
def load_solute(filename):
    if len(filename) < 3 :
        raise ValueError("Cannot determine filetype of input file '%s'" % filename)
    ext = filename[-3:]
    if ext=='mae' :
        molid = molecule.load('mae', filename)
    elif ext=='pdb' : # Need to convert to MAE because these things are stupid
        temp_mae = tempfile.mkstemp(suffix='.mae', prefix='dabble_input')[1]
        molid = molecule.load('pdb', filename)
        atomsel('all').write('mae', temp_mae)
        molecule.delete(molid)
        molid = molecule.load('mae', temp_mae)
    else :
        raise ValueError("Filetype '%s' currently unsupported for input protein" % ext) 
    return molid


# Accepts either a filename or a molecule id to get the dimensions of
def get_system_dimensions(molid=None, filename=None):
    if molid is not None :
        p = molecule.get_periodic(molid=molid)
    elif filename is not None :
        top = molecule.get_top()
        tmp_top = molecule.read(-1, 'mae', filename)
        p = molecule.get_periodic()
        molecule.delete(molecule.get_top())
        if top != -1: molecule.set_top(top)
    else:
        raise Exception('Specify molid or filename to get_system_dimensions')
    if p['a']==0.0 and p['b']==0.0 and p['c']==0.0:
        raise Exception('No periodic box found in membrane!')
    return p['a'], p['b'], p['c']


def concatenate_mae_files(output_filename, input_filenames=None, input_ids=None):
    """
    Concatenates several mae files together into one. Since this file format
    allows concatenation, the new file is a combined system of all
    the input molecules superimposed.
    Either takes a list of files to concatenate, or a list of VMD molecule
    ids. If molecule ids are given, VMD's interface is used to get the filename
    corresponding to each id.

    Args:
      output_filename (str): Filename to write
      input_filenames (list of str): List of input mae files to combine, OR:
      input_ids (list of int): List of input molecules to combine

    Returns:
      True if successful

    Raises:
      AssertionError : if there are no input files
    """

    if input_ids is not None:
      input_filenames = [ (molecule.get_filenames(id))[0] for id in input_ids ]
    assert len(input_filenames) > 0, 'need at least one input filename'
    outfile = open(output_filename, 'w')
    for line in open(input_filenames[0]):
        outfile.write(line)
    for input_filename in input_filenames[1:]:
        infile = open(input_filename)
        for i in xrange(5):
            infile.readline()
        for line in infile:
            outfile.write(line)
    outfile.close()
    return

# Combines molecules with given input ids
# Closes them, and returns the molecule id of a new molecule that combines them
# Puts that new molecule on top
def combine_molecules(input_ids) :
    output_filename = tempfile.mkstemp(suffix='.mae', prefix='dabble_combine')[1]
    concatenate_mae_files(output_filename, input_ids=input_ids);
    output_id = molecule.load('mae', output_filename)
    molecule.set_top(output_id)
    for i in input_ids: molecule.delete(i)
    return output_id

def write_ct_blocks(sel, output_filename):
    """
    Writes a mae format file containing the specified selection.

    Args:
      sel (str): the selection to write
      output_filename (str): the file to write to, including .mae extension

    Returns :
      length (int): the number of CT blocks written
    """
    users = sorted(set(atomsel(sel).get('user')))
    filenames = [(tempfile.mkstemp(suffix='.mae', prefix='dabble_tmp_user'))[1] for id in users]
    length = len(users)

    for id, fn in zip(users, filenames):
        tempsel = atomsel('user %f and (%s)' % (id, sel))
        sel2 = atomsel('index ' + ' '.join(map(str,set(tempsel.get('index')))))
        sel2.set('user', 0.0)
        sel2.write('mae', fn)

    # Option lets us specify if we should write a pdb/psf or just a mae file
    # Either way it writes a temp mae file, hacky but it works
    concatenate_mae_files(output_filename, input_filenames=filenames)

    # Clean up
    for filename in filenames:
        os.remove(filename) # delete temporary files
    return length

def tile_system(input_id, times_x, times_y, times_z):
# Read in the equilibrated bilayer file and put it as the active molceule
    new_resid = array(atomsel('all',molid=input_id).get('residue'))
    num_residues = new_resid.max()
    atomsel('all',molid=input_id).set('user', 2.)
    num_id_blocks = int(max(atomsel('all',molid=input_id).get('user')))
    wx, wy, wz = get_system_dimensions(molid=input_id)

# Move the lipids over, save that file, move them back, repeat, then stack all of those
# together to make a tiled membrane. Uses temporary mae files to save each "tile" since this
# format is easy to combine. Renumbers residues as it goes along.
    tile_filenames = []
    for nx in range(times_x):
        for ny in range(times_y):
            for nz in range(times_z):
                tx = array([nx * wx, ny * wy, nz * wz])
                atomsel('all',molid=input_id).moveby(tuple(tx))
                atomsel('all',molid=input_id).set('resid', new_resid)
                new_resid += num_residues
                (h,tile_filename) = tempfile.mkstemp(suffix='.mae', prefix='dabble_tile_tmp')
                tile_filenames.append(tile_filename)
                atomsel('all',molid=input_id).write('mae', tile_filename)
                atomsel('all',molid=input_id).moveby(tuple(-tx))

# Write all of these tiles together into one large bilayer
    (h,merge_output_filename) = tempfile.mkstemp(suffix='.mae', prefix='dabble_merge_tile_tmp')
    concatenate_mae_files(merge_output_filename, input_filenames=tile_filenames)

# Read that large bilayer file in as a new molecule and write it as the output file
    output_id = molecule.load('mae', merge_output_filename)
    molecule.set_periodic(output_id, -1, times_x * wx, times_y * wy, times_z * wz, 90.0, 90.0, 90.0)

# Rewrite user attribute
    for id in range(1, 1 + num_id_blocks * times_x * times_y * times_z) :
        atomsel('user %f' % id, molid=output_id).set('user', ((id-1) % num_id_blocks) + 1)
    #write_ct_blocks('all', output_filename)
    atomsel('all', molid=output_id).write('mae', merge_output_filename)

    #molecule.delete(molecule.get_top())
    #if top != -1: molecule.set_top(top)
    for tile_filename in tile_filenames:
        os.remove(tile_filename)
    #os.remove(merge_output_filename)
    return output_id


############################################################
#       SYSTEM PREP STRUCTURE MANIPULATION FUNCTIONS       #
############################################################


def copy_file(input_filename, output_filename):
    f = open(output_filename, 'w')
    for line in open(input_filename) :
       f.write(line)
    f.close() 

        
def tile_membrane_patch(input_id, min_xy_size, min_z_size):
    sys_dimensions = array([min_xy_size, min_xy_size, min_z_size])
    mem_dimensions = array(get_system_dimensions(molid=input_id))
    times_x, times_y, times_z = [int(times) for times in ceil(sys_dimensions / mem_dimensions)]
    if times_z >=2 :
      print "WARNING: you will need to add more solvent in the Z direction!"
      times_z = 1
    #assert times_z < 2, 'dabble currently does not support tiling in the z dimension'
    # add support for tiling in the z direction?
    if times_x == 1 and times_y == 1 and times_z == 1:
        output_id = input_id
        #copy_file(input_filename, output_filename)
    else:
        output_id = tile_system(input_id, times_x, times_y, times_z)
    return output_id, times_x, times_y, times_z


def set_cell_to_square_prism(xy_size, z_size):
    molecule.set_periodic(-1, -1, xy_size, xy_size, z_size, 90.0, 90.0, 90.0)


# Called on the entire lipid system to move it, then saves moved positions
# and reloads molecule in case the file needs to be concatenated
def center_system(molid, center_z=False):
    x, y, z = atomsel('all', molid=molid).center()

    # Move system so center is at origin or just xy plane?
    if center_z is True:
        atomsel('all', molid=molid).moveby((-x, -y, -z))
    else:
        atomsel('all', molid=molid).moveby((-x, -y, 0))

    # Save and reload the solute to record atom positions
    temp_mae = tempfile.mkstemp(suffix='.mae', prefix='dabble_centered')[1]
    atomsel('all',molid=molid).write('mae',temp_mae)
    molecule.delete(molid)
    new_id = molecule.load('mae', temp_mae)
    return new_id


def init_atoms():
    atomsel('all').set('beta', 1)


def atomsel_remaining(sel = 'all'):
    return atomsel('beta 1 and (%s)' % sel)


def num_atoms_remaining(sel = 'all'):
    return len(atomsel_remaining(sel))

def num_waters_remaining(water_sel = 'water and element O'):
    return len(atomsel_remaining(water_sel))

def num_lipids_remaining(lipid_sel):
    from numpy import unique
    return unique(atomsel_remaining(lipid_sel).get('fragment')).size


def remove_atoms(sel):
    s = atomsel('beta 1 and (%s)' % sel)
    s.set('beta', 0)
    return len(s)


def remove_residues(sel):
#    print 'DEBUG: remove_residues(%s) (%d atoms)' % (sel, len(atomsel(sel)))
    return remove_atoms('same residue as (%s)' % sel) # residue is guaranteed unique across insertions

def trim_water(wat_buffer, solute_sel):
    zcoord = atomsel(solute_sel).get('z')
    remove_residues('(not (%s)) and noh and z > %f' % (solute_sel, max(zcoord) + wat_buffer))
    remove_residues('(not (%s)) and noh and z < %f' % (solute_sel, min(zcoord) - wat_buffer))

def remove_z_residues(z_size, solute_sel):
    return remove_residues('(not (%s)) and noh and abs(z) > %f' % (solute_sel, z_size / 2.0))

def remove_xy_residues(xy_size, solute_sel, lipid_sel):
    half_xy_size = xy_size / 2.0
    box_sel_str = 'abs(x) > %f or abs(y) > %f' % (half_xy_size, half_xy_size)
    suspicious_lipid_residues=list(set(atomsel('(%s) and (%s)' % (lipid_sel, box_sel_str)).get('residue')))
    bad_lipids=list()
    for i in suspicious_lipid_residues:
       lipid_center=atomsel('noh and residue ' + str(i)).center()
       if abs(lipid_center[0]) > half_xy_size or abs(lipid_center[1]) > half_xy_size:
          bad_lipids.append(i)
    lipid_headgroup_sel='residue ' + ' '.join(map(str,bad_lipids))
    removal_sel_str = '(%s) or not (%s)' % (lipid_headgroup_sel, lipid_sel)
    total = remove_residues('noh and (%s) and (%s) and not (%s)' % (box_sel_str, removal_sel_str, solute_sel)) 
    return total


def remove_overlapping_residues(solute_sel,
                                lipid_sel,
                                lipid_friendly_sel=False,
                                dist=1.75,
                                lipid_dist=1.0, 
                                zh_mem_hyd=__MEMBRANE_HYDROPHOBIC_THICKNESS/2.0):
    clashing_sel = 'not (%s) and noh and not (%s) and (pbwithin %f of (noh and (%s)))' % (lipid_sel, solute_sel, dist, solute_sel)
    total = remove_residues(clashing_sel)
#    clashing_sel2 = '(%s) and noh and not (%s) and (pbwithin %f of (backbone and (%s)))' % (lipid_sel, solute_sel, dist, solute_sel)
#    total += remove_residues(clashing_sel2)
#    total += remove_residues('(%s) and not (%s)' % (clashing_sel, lipid_sel))
    if lipid_friendly_sel:
       clashing_sel_lipid = '(%s) and noh and not (%s) and pbwithin %f of (noh and (%s) and not (%s))' % (lipid_sel,solute_sel, dist, solute_sel,lipid_friendly_sel)
    else:
       clashing_sel_lipid = '(%s) and noh and not (%s) and pbwithin %f of (noh and (%s))' % (lipid_sel,solute_sel, dist, solute_sel)
    total += remove_residues(clashing_sel_lipid)
    return total


def remove_lipids_near_rings(solute_sel,
                             lipid_sel,
                             ring_sel='noh and resname HID HIE HIP HIS PHE TRP TYR and not backbone',
                             dist=1.75):
    total = 0
    solute_ring_sel = '(%s) and (%s)' % (solute_sel, ring_sel)
    total = remove_residues('noh and (%s) and not (%s) and pbwithin %f of (noh and (%s))' % (lipid_sel, solute_sel, dist, solute_ring_sel))
    return total

def remove_lipid_boundary_clash(pointy_type,ring_type,dist=1.0):
    total = 0
    sel1='noh and ' + pointy_type + ' and pbwithin ' + str(dist) + ' of noh and ' + ring_type
    sel2='noh and ' + ring_type   + ' and pbwithin ' + str(dist) + ' of noh and ' + pointy_type
    total = remove_residues('(' + sel1 + ')' + ' or ' + '(' + sel2 + ')')
    return total

def find_convertible_water_molecule(water_sel = 'water', min_ion_dist = 5.0):
    inclusion_sel = 'beta 1 and noh and (%s)' % water_sel
    exclusion_sel = 'beta 1 and not (%s)' % water_sel
    sel = atomsel('(%s) and not pbwithin %f of (%s)' % (inclusion_sel, min_ion_dist, exclusion_sel))
    assert len(sel) > 0, 'no convertible water molecules found in the selection "%s"' % str(sel)
    return sel.get('index')[random.randint(0, len(sel))]


def set_ion(gid, element):
    sel = atomsel('index %d' % gid)
    resname  = dict(Na='SOD', K='POT', Cl='CLA')[element]
    name     = dict(Na='NA', K='K', Cl='CL')[element]
    type     = dict(Na='NA', K='K', Cl='CL')[element]
    charge   = dict(Na=1, K=1, Cl=-1)[element]
    sel.set('element', element)
    sel.set('name', name)
    sel.set('type', type)
    sel.set('resname', resname)
    sel.set('chain', 'N')
    sel.set('segid', 'ION')
    sel.set('charge', charge)


def set_cations(element, filter_sel='none'):
    assert element in ['Na', 'K'], 'Supported cations are "Na" and "K"'
    for gid in tuple(atomsel('element K Na and not (%s)' % filter_sel)):
        set_ion(gid, element)


def convert_water_molecule_to_ion(gid, element):
    assert element in ['Na', 'K', 'Cl'], 'element must be "Na", "K", or "Cl"'
    element_received = atomsel('index %d' % gid).get('element')[0]
    assert element_received == 'O', 'Specify the oxygen atom of water molecule to convert (got %s)' % element_received
    remove_atoms('element H and same residue as index %d' % gid)
    set_ion(gid, element)
    return gid
    

def add_salt_ion(element):
    """
    Changes a water molecule to a salt ion.

    Args:
      element (str) : ion to add

    Raises:
      AssertionException : if element is not supported

    Returns:
      (int) the index of the water molecule that was replaced
    """
    assert element in ['Na', 'K', 'Cl'], 'element must be "Na", "K", or "Cl"'
    gid = find_convertible_water_molecule()
    convert_water_molecule_to_ion(gid, element)
    return gid


def write_final_system(opts, out_fmt, molid):
    """
    Writes the final output in whatever format(s) are requested.
    Always writes a mae format file

    Args:
      opts (argparse) : options passed to dabble
      out_fmt (str): format to write the output to
      molid (int): VMD molecule_id to write

    Returns:
      (str) main final filename written
    """

    # Write a mae file always, removing the prefix from the output file
    mae_name = '.'.join(map(str,opts.output_filename.rsplit('.')[:-1])) + '.mae'
    write_ct_blocks(sel='beta 1', output_filename=mae_name)

    # If a converted output format (pdb or dms) desired, write that here
    # and the mae is a temp file that can be deleted
    if out_fmt=='dms' :
        temp_mol = molecule.load('mae',mae_name)
        atomsel('all',molid=temp_mol).write(out_fmt, opts.output_filename)
        molecule.delete(temp_mol)
        os.remove(mae_name)

    # For pdb, write an AMBER leap compatible pdb, don't trust the VMD
    # pdb writing routine
    if out_fmt=='pdb' : 
        #import dabbleparam
        temp_mol = molecule.load('mae', mae_name)
        atomsel('all',molid=temp_mol).write(out_fmt, opts.out_filename)
        #dabbleparam.write_amber_pdb(opts.output_filename, molid=temp_mol)
        molecule.delete(temp_mol)

    # If we want a parameterized format like amber or charmm, a psf must
    # first be written which does the atom typing, etc
    if out_fmt=='charmm' or out_fmt=='amber' :
        import dabbleparam
        temp_mol = molecule.load('mae', mae_name)
        write_psf_name = mae_name.replace('.mae','')
        topos = dabbleparam.write_psf(write_psf_name, molid=temp_mol, lipid_sel=opts.lipid_sel)

    # For amber format files, invoke the parmed chamber routine
    if out_fmt=='amber' :
        print("\nINFO: Writing AMBER format files with CHARMM parameters. This may take a moment...\n")
        dabbleparam.psf_to_amber(write_psf_name, write_psf_name, topos, molid) 

    return opts.output_filename


def check_write_ok(filename, out_fmt, overwrite=False):
    """
    Checks if the output files for the requested format exists,
    and prints out an error message if the current options
    don't allow overwriting them.

    Args:
      filename (str) : Output filename requested
      out_fmt (str) : Output format requested. All intermediate
      files involved in writing to this format will be checked for
      existence.
      overwrite (bool) : True if overwriting is allowed

    Returns:
      True if it okay to overwrite
      Quits the program otherwise
    """
    if overwrite is True :
        return True
    
    # Generate file suffixes to search for
    prefix = '.'.join(filename.split('.')[:-1])
    suffixes = ['mae']
    if out_fmt=='dms':
        suffixes.append('dms')
    elif out_fmt=='pdb':
        suffixes.append('pdb')
    elif out_fmt=='charmm':
        suffixes.extend(['psf','pdb'])
    elif out_fmt=='amber':
        suffixes.extend(['psf','pdb','prmtop','inpcrd'])
            
    exists = []
    for s in suffixes:
        if os.path.isfile('%s.%s'%(prefix,s)):
            exists.append('%s.%s'%(prefix,s))

    if len(exists):
        print("\nERROR: The following files exist and would be overwritten:\n")
        print(  "       %s\n" % ' '.join(exists))
        print(  "       Won't overwrite, exiting.")
        print(  "       Run with -O to overwrite files next time.")
        quit(1)

    return False
