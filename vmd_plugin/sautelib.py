# sautelib.py
# 
# $Date: 2012/01/27 $
# Dan Arlow and Thomas Mildorf
# Dan.Arlow@deshawresearch.com
# Thomas.Mildorf@deshawresearch.com
# D E Shaw Research LLC



import os
import tempfile

from numpy import *

import vmd
from atomsel import atomsel

# Yes, these are different modules. With different functions. And we need both.
# VMD is the worst.
from Molecule import * # http://www.ks.uiuc.edu/Research/vmd/current/ug/node188.html
from molecule import * # http://www.ks.uiuc.edu/Research/vmd/current/ug/node182.html


############################################################
#                       CONSTANTS                          #
############################################################


__MEMBRANE_HYDROPHOBIC_THICKNESS = 30.0
__MEMBRANE_FULL_THICKNESS = 50.0

__1M_SALT_IONS_PER_WATER = 0.018

random.seed(2011)

############################################################
#             STRUCTURE ANALYSIS FUNCTIONS                 #
############################################################


# default sel was solute_sel = 'ctnumber 1'
# seemed like an abstraction violation
# 
def get_net_charge(sel):
    charge = array(atomsel(sel).get('charge'))
    if charge.size == 0:
        return 0
    assert any(charge != 0), 'all charges are zero! did you use a maeff file?'
    net_charge = sum(charge)
    rslt = round(net_charge)
    assert abs(rslt - net_charge) < .01, 'total charge is not integral (within a tolerance of .01) -- something wrong with your force field?'
    return int(rslt)

def get_system_net_charge():
    return get_net_charge('beta 1')


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


# former default:  solute_sel = 'ctnumber 1'
def solute_xy_diameter(solute_sel):
    a = atomsel(solute_sel)
    return diameter(transpose([a.get('x'), a.get('y')]))

# Accepts either a molecule id or file to load. 
def get_solute_sel(molid=None, filename=None):
    atomsel('all').set('user',1.)
    if molid is not None:
        #sel = 'ctnumber ' + ' '.join(map(str, set(atomsel('all', molid=molid).get('ctnumber'))))
        sel = 'serial ' + ' '.join(map(str, set(atomsel('all', molid=molid).get('serial'))))
    elif filename is not None:
        top = molecule.get_top()
        tmp_top = molecule.read(-1, 'mae', input_filename)
        sel = 'serial ' + ' '.join(map(str, set(atomsel('all').get('serial'))))
        molecule.delete(tmp_top)
        if top != -1: molecule.set_top(top)
    else:
      raise Exception('Specify molid or filename to get_solute_sel')

    # Set the user attribute to 1 in case ctnumber isn't an available attribute
    return sel


def move(solute_sel, vector, molid=-1):
    atomsel(solute_sel, molid=molid).moveby(vector)


# former default:  solute_sel = 'ctnumber 1'
def get_cell_size(mem_buf,
                  wat_buf,
                  solute_sel,
                  molid=None,
                  filename=None,
                  zh_mem_full=__MEMBRANE_FULL_THICKNESS / 2.0,
                  zh_mem_hyd=__MEMBRANE_HYDROPHOBIC_THICKNESS / 2.0,
                  z_move=0):

    if filename is not None :
        top = molecule.get_top()    
        tmp_top = molecule.read(-1, 'mae', filename)

    if molid is None:
        molid = molecule.get_top()

    move(solute_sel, (0, 0, z_move), molid=molid)

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
    try:
        abs(get_net_charge(str(cations))-num_cations) > 0.01
    except:
        # Check for bonded cations
        nonbonded_cation_index=[atomsel('index %d' %x).get('index')[0] for x in nonzero(array(map(len,atomsel('element %s' %cation).bonds)) == 0)[0]]
        if len(nonbonded_cation_index)==0:
            cations=atomsel('none')
        else:
            cations=atomsel_remaining('index %s' %' '.join(map(str,nonbonded_cation_index)))
        num_cations = len(cations)
        assert abs(get_net_charge(str(cations))-num_cations) < 0.01, 'num cations and net cation charge are not equal'
    try:
        abs(get_net_charge(str(anions))+num_anions) > 0.01
    except:
        # Check for bonded anions
        non_bonded_anion_index=[atomsel('index %d' %x).get('index')[0] for x in nonzero(array(map(len,atomsel_remaining('element %s' %anion).bonds)) == 0)[0]]
        if len(non_bonded_anion_index)==0:
            anions=atomsel('none')
        else:
            anions=atomsel_remaining('index %s' %' '.join(map(str,nonbonded_anion_index)))
        num_anions = len(anions)
        assert abs(get_net_charge(str(anions))+num_anions) < 0.01, 'num anions and abs anion charge are not equal'
    num_waters = num_atoms_remaining(water_sel)
    num_for_conc = int(round(__1M_SALT_IONS_PER_WATER * num_waters * conc))
    pos_ions_needed = num_for_conc - num_cations
    neg_ions_needed = num_for_conc - num_anions
    system_charge = get_system_net_charge()

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
    return p['a'], p['b'], p['c']


def concatenate_mae_files(output_filename, input_filenames):
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


def read_combined_mae_file(input_filenames):
    assert len(input_filenames) > 0, 'need at least one input filename'
    tmp_filename = tempfile.mktemp(suffix='mae', prefix='saute_tmp')
    concatenate_mae_files(tmp_filename, input_filenames)
    m = molecule.read(-1, 'mae', tmp_filename)
    os.remove(tmp_filename)
    return m


def write_ct_blocks(sel, output_filename):
    # Get all ctnumbers or use user field 
    try :
        ctnumbers = sorted(set(atomsel(sel).get('ctnumber')))
        filenames = [tempfile.mktemp(suffix='mae', prefix='saute_tmp') for ctnumber in ctnumbers]
        length = len(ctnumbers)
        
        for ctnumber, filename in zip(ctnumbers, filenames):
            atomsel('ctnumber %d and (%s)' % (ctnumber, sel)).write('mae', filename)
        concatenate_mae_files(output_filename, filenames)
    except ValueError :
        users = sorted(set(atomsel(sel).get('user')))
        filenames = [tempfile.mktemp(suffix='mae', prefix='saute_tmp_user') for id in users]
        length = len(users)

        for id, fn in zip(users, filenames):
            tempsel = atomsel('user %f and (%s)' % (id, sel))
            sel2 = atomsel('index ' + ' '.join(map(str,set(tempsel.get('index')))))
            sel2.set('user', 0.0)
            sel2.write('mae', fn)

    concatenate_mae_files(output_filename, filenames)
    for filename in filenames:
        os.remove(filename) # delete temporary files
    return length

def tile_system(input_filename, output_filename, times_x, times_y, times_z):
# Read in the equilibrated bilayer file and put it as the active molceule
    top = molecule.get_top()
    tmp_top = molecule.read(-1, 'mae', input_filename)
    new_resid = array(atomsel('all').get('residue'))
    num_residues = new_resid.max()
    atomsel('all').set('user', 2.)
    try :
        print("Found and using CT number")
        using_ct = True
        num_ct_blocks = max(atomsel('all').get('ctnumber')) # specific to mae files
    except ValueError : # no ctnumber attribute
        num_id_blocks = int(max(atomsel('all').get('user')))
        using_ct = False
    wx, wy, wz = get_system_dimensions(filename=input_filename)

# Move the lipids over, save that file, move them back, repeat, then stack all of those
# together to make a tiled membrane. Uses temporary mae files to save each "tile" since this
# format is easy to combine. Renumbers residues as it goes along.
    tile_filenames = []
    for nx in range(times_x):
        for ny in range(times_y):
            for nz in range(times_z):
                tx = array([nx * wx, ny * wy, nz * wz])
                atomsel('all').moveby(tuple(tx))
                atomsel('all').set('resid', new_resid)
                new_resid += num_residues
                tile_filename = tempfile.mktemp(suffix='mae', prefix='saute_tile_tmp')
                tile_filenames.append(tile_filename)
                atomsel('all').write('mae', tile_filename)
                atomsel('all').moveby(tuple(-tx))
    molecule.delete(tmp_top)

# Write all of these tiles together into one large bilayer
    merge_output_filename = tempfile.mktemp(suffix='mae', prefix='saute_merge_tile_tmp')
    concatenate_mae_files(merge_output_filename, tile_filenames)

# Read that large bilayer file in as a new molecule and write it as the output file
    molecule.read(-1, 'mae', merge_output_filename)
    molecule.set_periodic(-1, -1, times_x * wx, times_y * wy, times_z * wz, 90.0, 90.0, 90.0)

# Rewrite ctnumber. Will not run loop if ctnumber isn't present in this file type
    if using_ct :
        for ctnumber in range(1, 1 + num_ct_blocks * times_x * times_y * times_z) :
            atomsel('ctnumber %d' % ctnumber).set('ctnumber', ((ctnumber-1) % num_ct_blocks) + 1)
    else :
        for id in range(1, 1 + num_id_blocks * times_x * times_y * times_z) :
            atomsel('user %f' % id).set('user', ((id-1) % num_id_blocks) + 1)
    write_ct_blocks('all', output_filename)

    molecule.delete(molecule.get_top())
    if top != -1: molecule.set_top(top)
    for tile_filename in tile_filenames:
        os.remove(tile_filename)
    os.remove(merge_output_filename)
    return


############################################################
#       SYSTEM PREP STRUCTURE MANIPULATION FUNCTIONS       #
############################################################


def copy_file(input_filename, output_filename):
    f = open(output_filename, 'w')
    for line in open(input_filename) :
       f.write(line)
    f.close() 

        
def tile_membrane_patch(input_filename, output_filename, min_xy_size, min_z_size):
    sys_dimensions = array([min_xy_size, min_xy_size, min_z_size])
    mem_dimensions = array(get_system_dimensions(filename=input_filename))
    times_x, times_y, times_z = [int(times) for times in ceil(sys_dimensions / mem_dimensions)]
    assert times_z < 2, 'saute currently does not support tiling in the z dimension'
    # add support for tiling in the z direction?
    if times_x == 1 and times_y == 1 and times_z == 1:
        copy_file(input_molid, output_filename)
    else:
        tile_system(input_filename, output_filename, times_x, times_y, times_z)
    return times_x, times_y, times_z


def set_cell_to_square_prism(xy_size, z_size):
    molecule.set_periodic(-1, -1, xy_size, xy_size, z_size, 90.0, 90.0, 90.0)


# former default: lipid_system_sel = 'not ctnumber 1'
def center_membrane_system(solute_sel, lipid_sel, new_z_center=0):
    #    lipid_system_sel = 'not ctnumber 1'
    lipid_system_sel = 'not (%s)' % str(solute_sel)
    x, y, z = atomsel('%s' % lipid_sel).center()
    atomsel(lipid_system_sel).moveby((-x, -y, new_z_center - z))
    return x, y, z


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


# former default: solute_sel = 'not ctnumber 1'
def remove_z_residues(z_size, solute_sel):
    return remove_residues('(not (%s)) and noh and abs(z) > %f' % (solute_sel, z_size / 2.0))

# former default: solute_sel = 'not ctnumber 1'
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
#    ctnumber = dict(Na=4, K=4, Cl=5)[element]
    charge   = dict(Na=1, K=1, Cl=-1)[element]
    sel.set('element', element)
    sel.set('name', name)
    sel.set('type', type)
    sel.set('resname', resname)
    sel.set('chain', 'N')
    sel.set('segid', 'ION')
#    sel.set('ctnumber', ctnumber)
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
    assert element in ['Na', 'K', 'Cl'], 'element must be "Na", "K", or "Cl"'
    gid = find_convertible_water_molecule()
    convert_water_molecule_to_ion(gid, element)
    return gid


def write_remaining_atoms(output_filename):
    write_ct_blocks('beta 1', output_filename)
    return num_atoms_remaining()
    

