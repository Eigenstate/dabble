# This does NOT read any parameters from the file, just the
# coordinates and velocities.
# 12-16-2012 Markus Dittrich: Removed code that auto generates the
#            vdw_combined section since Anton treats these as nbfixes
#            and subsequently fails due to their large number

import optparse
import math
import sys
import re
import copy
from schrodinger import structure, structureutil
from schrodinger.infra import mm
import schrodinger.structutils.assignbondorders as assign_bo
import numpy, re

blocks = {}
amber_st = None

def createMassTable():
    mass_table = {}
    mass_table[math.floor(0.5 + 100*0.0000)] = 0          #DU
    mass_table[math.floor(0.5 + 100*1.0079)] = 1          #H 
    mass_table[math.floor(0.5 + 100*2.0136)] = 1          #D
    mass_table[math.floor(0.5 + 100*4.0026)] = 2          #He
    mass_table[math.floor(0.5 + 100*6.941)] = 3           #Li
    mass_table[math.floor(0.5 + 100*9.0122)] = 4          #Be
    mass_table[math.floor(0.5 + 100*10.811)] = 5          #B 
    mass_table[math.floor(0.5 + 100*12.011)]  = 6         #C 
    mass_table[math.floor(0.5 + 100*14.0067)] = 7         #N 
    mass_table[math.floor(0.5 + 100*15.9994)] = 8         #O 
    mass_table[math.floor(0.5 + 100*18.9984)] = 9         #F 
    mass_table[math.floor(0.5 + 100*20.1797)] = 10        #Ne
    mass_table[math.floor(0.5 + 100*22.9898)] = 11        #Na
    mass_table[math.floor(0.5 + 100*24.3050)] = 12        #Mg
    mass_table[math.floor(0.5 + 100*26.9815)] = 13        #Al
    mass_table[math.floor(0.5 + 100*28.0855)] = 14        #Si
    mass_table[math.floor(0.5 + 100*30.9738)] = 15        #P 
    mass_table[math.floor(0.5 + 100*32.064)] = 16         #S 
    mass_table[math.floor(0.5 + 100*35.4527)] = 17        #Cl
    mass_table[math.floor(0.5 + 100*39.948)] = 18         #Ar
    mass_table[math.floor(0.5 + 100*39.0983)] = 19        #K 
    mass_table[math.floor(0.5 + 100*40.078)] = 20         #Ca
    mass_table[math.floor(0.5 + 100*44.955912)] = 21      #Sc
    mass_table[math.floor(0.5 + 100*47.867)] = 22         #Ti
    mass_table[math.floor(0.5 + 100*50.9415)] = 23        #V 
    mass_table[math.floor(0.5 + 100*51.9961)] = 24        #Cr
    mass_table[math.floor(0.5 + 100*54.938045)] = 25      #Mn
    mass_table[math.floor(0.5 + 100*55.845)] = 26         #Fe
    mass_table[math.floor(0.5 + 100*58.933195)] = 27      #Co
    mass_table[math.floor(0.5 + 100*58.6934)] = 28        #Ni
    mass_table[math.floor(0.5 + 100*63.546)] = 29         #Cu
    mass_table[math.floor(0.5 + 100*65.409)] = 30         #Zn
    mass_table[math.floor(0.5 + 100*69.723)] = 31         #Ga
    mass_table[math.floor(0.5 + 100*72.64)] = 32          #Ge
    mass_table[math.floor(0.5 + 100*74.92160)] = 33       #As
    mass_table[math.floor(0.5 + 100*78.96)] = 34          #Se
    mass_table[math.floor(0.5 + 100*79.904)] = 35         #Br
    mass_table[math.floor(0.5 + 100*85.4678)] = 37        #Rb
    mass_table[math.floor(0.5 + 100*87.62)] = 38          #Sr
    mass_table[math.floor(0.5 + 100*88.90585)] = 39       #Y 
    mass_table[math.floor(0.5 + 100*91.224)] = 40         #Zr
    mass_table[math.floor(0.5 + 100*92.90638)] = 41       #Nb
    mass_table[math.floor(0.5 + 100*95.96)] = 42          #Mo
    mass_table[math.floor(0.5 + 100*97.9072)] = 43        #Tc
    mass_table[math.floor(0.5 + 100*101.07)] = 44         #Ru
    mass_table[math.floor(0.5 + 100*102.90550)] = 45      #Rh
    mass_table[math.floor(0.5 + 100*106.42)] = 46         #Pd
    mass_table[math.floor(0.5 + 100*107.8682)] = 47       #Ag
    mass_table[math.floor(0.5 + 100*112.411)] = 48        #Cd
    mass_table[math.floor(0.5 + 100*118.710)] =  49       #Sn
    mass_table[math.floor(0.5 + 100*121.75)] = 51         #Sb
    mass_table[math.floor(0.5 + 100*126.90447)] = 53      #I 
    mass_table[math.floor(0.5 + 100*131.30)] = 54         #Xe
    mass_table[math.floor(0.5 + 100*132.9054519)] = 55    #Cs
    mass_table[math.floor(0.5 + 100*150.35)] = 62         #Sm 
    mass_table[math.floor(0.5 + 100*151.25)] = 63         #Eu 
    mass_table[math.floor(0.5 + 100*183.84)] = 74         #W 
    mass_table[math.floor(0.5 + 100*186.2)] = 75          #Re 
    mass_table[math.floor(0.5 + 100*190.2)] = 76          #Os
    mass_table[math.floor(0.5 + 100*195.084)] = 78        #Pt
    mass_table[math.floor(0.5 + 100*196.966569)] = 79     #Au
    mass_table[math.floor(0.5 + 100*200.59)] = 80         #Hg
    mass_table[math.floor(0.5 + 100*207.2)] = 82          #Pb
    mass_table[math.floor(0.5 + 100*208.98040)] = 83      #Bi
    mass_table[math.floor(0.5 + 100*231.0)] = 92          #U

    return mass_table

def parsePrmtop(ifname):
    f = open(ifname)
    lines = f.readlines()

    blks = {}
    
    re_flag = re.compile(r'\%FLAG\s+(\S+)', re.IGNORECASE)
    re_format = re.compile(r'\%FORMAT\s*\((\d+)\w(\d+)\S*\)', re.IGNORECASE)
    re_comment = re.compile(r'\%\S+', re.IGNORECASE)
    flag_name = ''
    flag_length = 0
    for l in lines:
        flag_match = re_flag.search(l)
        format_match = re_format.search(l)
        comment_match = re_comment.search(l)
        if flag_match:
            flag_name = flag_match.group(1)
            flag_size = 0
            flag_length = 0
            blks[flag_name] = []
        elif format_match:
            flag_size = int(format_match.group(1))
            flag_length = int(format_match.group(2))
        elif comment_match:
            continue
        else:
            if not (flag_name and flag_size and flag_length):
                print 'cannot recognize flag.'
            for i in range(flag_size):
                start = i*flag_length
                end = start + flag_length
                if end > len(l):
                    end = len(l)
                element = l[start:end].strip()
                if element:
                    blks[flag_name].append(element)
                else:
                    break
    return blks

def convertTop2Ffio(ofname):
    global blocks

    charge = blocks['CHARGE']
    mass = blocks['MASS']
    type = blocks['AMBER_ATOM_TYPE']
    type_index = blocks['ATOM_TYPE_INDEX']
    ntype = len(set(type_index))
    vdw_type_symbol = range(ntype)
    for t, ti in zip(type, type_index):
        vdw_type_symbol[int(ti)-1] = t

    s = """
  ffio_ff {
    s_ffio_name
    s_ffio_comb_rule
    i_ffio_version
    :::
    AMBER
    ARITHMETIC/GEOMETRIC
    1"""

    s += """
    ffio_sites[%d] { 
      s_ffio_type
      r_ffio_charge
      r_ffio_mass
      s_ffio_vdwtype
      :::\n""" % len(charge)
    i = 1
    for (c, m) in zip(charge, mass):
        ti = type_index[i-1]
        t = vdw_type_symbol[int(ti)-1]
        s +=  '      %d  atom  %f   %s   %s\n' % ( i, float(c)/18.2223, m, t)
        i += 1
    s += '      :::\n'
    s += '    }'         

    
    bond = copy.deepcopy(blocks['BONDS_WITHOUT_HYDROGEN'])
    bond.extend(blocks['BONDS_INC_HYDROGEN'])
    row = len(bond) / 3
    bond = numpy.array(bond)
    bond = bond.reshape(row, 3)
    bond_force = blocks['BOND_FORCE_CONSTANT']
    bond_dist = blocks['BOND_EQUIL_VALUE']

    s += """
    ffio_bonds[%d] { 
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::\n""" % row
    i = 1
    for a in bond:
        ai = int(a[0]) / 3 + 1
        aj = int(a[1]) / 3 + 1
        index = int(a[2]) - 1
        r0 = bond_dist[index]
        k = bond_force[index]
        s += '      %d  %d  %d   Harm   %s  %s\n' % ( i, ai, aj, r0, k)
        i += 1
    s += '      :::\n'
    s +=  '    }'         

    angle = copy.deepcopy(blocks['ANGLES_WITHOUT_HYDROGEN'])
    angle.extend(blocks['ANGLES_INC_HYDROGEN'])
    row = len(angle) / 4
    angle = numpy.array(angle)
    angle = angle.reshape(row, 4)
    angle_force = blocks['ANGLE_FORCE_CONSTANT']
    angle_value = blocks['ANGLE_EQUIL_VALUE']

    s += """
    ffio_angles[%d] { 
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::\n""" % len(angle)
    i = 1
    for a in angle:
        ai = int(a[0]) / 3 + 1
        aj = int(a[1]) / 3 + 1
        ak = int(a[2]) / 3 + 1
        index = int(a[3]) - 1
        r0 = math.degrees(float(angle_value[index]))
        k = angle_force[index]
        s += '      %d  %d  %d  %d  Harm   %f  %s\n' % ( i, ai, aj, ak, r0, k)
        i += 1
    s += '      :::\n'
    s += '    }'
    
    dihedral = copy.deepcopy(blocks['DIHEDRALS_WITHOUT_HYDROGEN'])
    dihedral.extend(blocks['DIHEDRALS_INC_HYDROGEN'])
    row = len(dihedral) / 5
    dihedral = numpy.array(dihedral)
    dihedral = dihedral.reshape(row, 5)
    dihedral_force = blocks['DIHEDRAL_FORCE_CONSTANT']
    dihedral_period = blocks['DIHEDRAL_PERIODICITY']
    dihedral_phase = blocks['DIHEDRAL_PHASE']

    s += """
    ffio_dihedrals[%d] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      i_ffio_al
      s_ffio_funct
      r_ffio_c0
      r_ffio_c1
      r_ffio_c2
      r_ffio_c3
      r_ffio_c4
      r_ffio_c5
      r_ffio_c6
      r_ffio_c7
      :::\n""" % len(dihedral)

    i = 1
    for a in dihedral:
        ai = math.fabs(int(a[0])) / 3 + 1
        aj = math.fabs(int(a[1])) / 3 + 1
        ak = math.fabs(int(a[2])) / 3 + 1
        al = math.fabs(int(a[3])) / 3 + 1
        index = int(a[4]) - 1
        k = dihedral_force[index]
        n = int(float(dihedral_period[index]))
        p = math.degrees(float(dihedral_phase[index]))
        k = dihedral_force[index]
        funct = 'Proper_Trig'
        if (int(a[3]) < 0):
           funct = 'Improper_Trig'
        s += '      %d  %d  %d  %d  %d %s %f %s ' % ( i, ai, aj, ak, al, funct, p, k)
        for j in range(1, 7):
            if j == n:
                s += '%s ' % k
            else:
                s += '0.0 '
        s += '\n'
        i += 1
    s += '      :::\n'
    s += '    }'

    pair = []
    for a in dihedral:
        ai = int(a[0])
        aj = int(a[1])
        ak = int(a[2])
        al = int(a[3])
        if ak < 0 or al < 0:
            continue
        pair.append((ai/3 + 1, al/3 + 1))

    s += """
    ffio_pairs[%d] {
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      :::\n""" % (len(pair)*2)

    i = 1
    for a in pair:
        s += '      %d  %d %d Coulomb 0.8333\n' % (i, a[0], a[1])
        i += 1
    for a in pair:
        s += '      %d  %d %d LJ 0.5\n' % (i, a[0], a[1])
        i += 1
    s += '      :::\n'
    s += '    }'

    num_exclusion = blocks['NUMBER_EXCLUDED_ATOMS']
    excluded_atom_list = blocks['EXCLUDED_ATOMS_LIST']
    exclusion = []
    k = 0
    for i in range(len(num_exclusion)):
        for j in range(int(num_exclusion[i])):
            if excluded_atom_list[k] != '0':
                exclusion.append((str(i+1), excluded_atom_list[k]))
            k += 1
                         
    s += """
    ffio_exclusions[%d] { 
      i_ffio_ai
      i_ffio_aj
      :::\n""" % len(exclusion)
    i = 1
    for a in exclusion:
        s += '      %d  %s  %s\n' % (i, a[0], a[1])
        i += 1
    s += '      :::\n'
    s += '    }'

    acoeff = blocks['LENNARD_JONES_ACOEF']
    bcoeff = blocks['LENNARD_JONES_BCOEF']
    nonbonded = blocks['NONBONDED_PARM_INDEX']

    vdwtype = []
    vdwtype_combined = []
    k = 0
    for i in range(ntype):
        for j in range(ntype):
            index = int(nonbonded[k]) - 1
            A = float(acoeff[index])
            B = float(bcoeff[index])
            sigma = 0.0
            epsilon = 0.0
            if A != 0 or B != 0:
                epsilon = B*B/4.0/A
                sigma = (A/B)**(1.0/6.0)
            #print i, j, k, type_list[i], type_list[j], acoeff[index], bcoeff[index]
            #vdwtype_combined.append((vdw_type_symbol[i], vdw_type_symbol[j], sigma, epsilon))
            if i == j:
                vdwtype.append((vdw_type_symbol[i], sigma, epsilon))
            k += 1
    
    s += """
    ffio_vdwtypes[%d] {
      s_ffio_name
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::\n""" % ntype

    i = 1
    for a in vdwtype:
        s += '      %d  %s  LJ12_6_sig_epsilon %f %f\n' % (i, a[0], a[1], a[2])
        i += 1
    s += '      :::\n'
    s += '    }\n'         

    #s += """
    #ffio_vdwtypes_combined[%d] {
    #  s_ffio_name1
    #  s_ffio_name2
    #  s_ffio_funct
    #  r_ffio_c1
    #  r_ffio_c2
    #  :::\n""" % (ntype*ntype)

    #i = 1
    #for a in vdwtype_combined:
    #    s += '      %d  %s  %s LJ12_6_sig_epsilon %f %f\n' % (i, a[0], a[1], a[2], a[3])
    #    i += 1
    #s += '      :::\n'
    #s += '    }\n'         

    return s

def convertCrd2Mae(ifname, ofname):
    global blocks, amber_st
    f = open(ifname)
    lines = f.readlines()
    title = lines[0]
    try:
        natom = int(lines[1].split()[0])
    except:
        print '%s is not compatible with the amber crd format.' % ifname
        print 'Exiting...'
        sys.exit(1)

    ct = mm.mmct_ct_new(natom)
    st = structure.Structure(ct)
    box_line = lines[-1]
    format = re.compile(r'(\d+)?[edf](\d+)\.(\d+)%', re.I)
    coor =  [ s[i:i+12].strip() for s in lines[2:-1] for i in range(0, len(s), 12) ]
    coor = filter(None, coor)
    # Sanity check
    if len(coor)/6 != natom :
        print("ERROR: Read coordinates and velocities for %f atoms, but have %d atoms" % (len(coor)/6.,natom))
        print("       Check the format of your restart file. It should be 6F12.7 on each line")
        quit(1)
    row = len(coor)/3
    coor = numpy.array(coor)
    coor = coor.reshape(row, 3)

    atom_name = blocks['ATOM_NAME']
    mass = blocks['MASS']
    mass_table = createMassTable()
    for i in range(natom):
        st.atom[i+1].pdbname = atom_name[i]
        st.atom[i+1].atom_name = atom_name[i]
        mass_index = math.floor(0.5 + 100*float(mass[i]))
        try:
            st.atom[i+1].atomic_number = mass_table[mass_index]
        except:
            print 'There is no atomic number for atom %d and its mass, %s.\n0 is assigned to this atom.' % (i+1, mass[i])
            st.atom[i+1].atomic_number = 0
            
        wild_type = mm.mmat_get_wildcard( st.atom[i+1].atomic_number)
        new_color = mm.mmat_get_color( wild_type )
        mm.mmct_atom_set_color( st.handle, i+1, new_color )

        st.atom[i+1].x = float(coor[i][0])
        st.atom[i+1].y = float(coor[i][1])
        st.atom[i+1].z = float(coor[i][2])
        st.atom[i+1].property['r_ffio_x_vel'] = float(coor[natom+i][0])*20.4550
        st.atom[i+1].property['r_ffio_y_vel'] = float(coor[natom+i][1])*20.4550
        st.atom[i+1].property['r_ffio_z_vel'] = float(coor[natom+i][2])*20.4550


    residue_name = blocks['RESIDUE_LABEL']
    residue_number = blocks['RESIDUE_POINTER']
    residue_number.append(natom+1)
    resnum = 0
    for i, a in enumerate(st.atom):
        if a.index >= int(residue_number[resnum]):
            resnum += 1
        a.resnum = resnum
        a.pdbres = residue_name[resnum - 1]

    bond = copy.deepcopy(blocks['BONDS_WITHOUT_HYDROGEN'])
    bond.extend(blocks['BONDS_INC_HYDROGEN'])
    row = len(bond) / 3
    bond = numpy.array(bond)
    bond = bond.reshape(row, 3)
    for a in bond:
        i = int(a[0])/3 + 1
        j = int(a[1])/3 + 1
        if st.atom[i].atomic_number > 1 or st.atom[j].atomic_number > 1:
            st.atom[i].addBond(j, 1)
    st.retype()
    print("Assigning bond order (will take a long time)")
    assign_bo.assign_st(st)

    ax = 100.0
    by = 100.0
    cz = 100.0
    box = box_line.split()
    if len(box) == 6:
        ax = float(box[0])
        by = float(box[1])
        cz = float(box[2])
    st.property['r_chorus_box_ax'] = ax
    st.property['r_chorus_box_ay'] = 0.0
    st.property['r_chorus_box_az'] = 0.0
    st.property['r_chorus_box_bx'] = 0.0
    st.property['r_chorus_box_by'] = by
    st.property['r_chorus_box_bz'] = 0.0
    st.property['r_chorus_box_cx'] = 0.0
    st.property['r_chorus_box_cy'] = 0.0
    st.property['r_chorus_box_cz'] = cz

    st.property['s_ffio_ct_type'] = 'solute'
    print("Writing system")
    st.write(ofname)
    amber_st = st
    return ''.join(open(ofname).readlines())

def buildConstraints():
    global blocks, amber_st
    st = amber_st

    natom = st.atom_total
    bond = blocks['BONDS_INC_HYDROGEN']
    row = len(bond) / 3
    bond = numpy.array(bond)
    bond = bond.reshape(row, 3)
    bond_dist = blocks['BOND_EQUIL_VALUE']

    bond_table = {}
    for a in bond:
        ai = int(a[0]) / 3 + 1
        aj = int(a[1]) / 3 + 1
        key1 = str(ai) + ' ' + str(aj)
        key2 = str(aj) + ' ' + str(ai)
        index = int(a[2]) - 1
        r0 = bond_dist[index]
        bond_table[key1] = r0
        bond_table[key2] = r0

    bond_constraints = []
    functs = []
    for a in st.atom:
        atom_constr = []
        if a.atomic_number > 1:
            for b in a.bonded_atoms:
                if b.atomic_number == 1:
                    atom_constr.append(b.index)
        bond_constraints.append(atom_constr)

    i = 1
    for a in bond_constraints:
        nconstr = len(a)
        if nconstr == 2 and st.atom[i].atomic_number == 8:
            functs.append('HOH')
        else:
            fct = 'AH%d' % nconstr
            functs.append(fct)
        i += 1
    constraints = []
    for (i, a) in zip(range(natom), bond_constraints):
        if functs[i] == 'AH0' or functs[i] == 'HOH':
            continue
        temp_constr = [ 0, 0, 0, 0, 0, '', 0.0, 0.0, 0.0, 0.0, 0.0]
        temp_constr[0] = i+1
        temp_constr[5] = functs[i]
        for j in range(len(a)):
            temp_constr[j+1] = a[j]
            key = str(i+1) + ' ' + str(a[j])
            temp_constr[j+6] = float(bond_table[key])
        constraints.append(temp_constr)

    for i in range(natom):
        if functs[i] == 'HOH':
            a = st.atom[i+1]
            hs = []
            for b in a.bonded_atoms:
                hs.append(b.index)
            if len(hs) != 2:
                print 'this is not water.'
                break
            key1 = str(a.index) + ' ' + str(hs[0])
            key2 = str(a.index) + ' ' + str(hs[1])
            key3 = str(hs[0]) + ' ' + str(hs[1])
            oh1 = float(bond_table[key1])
            oh2 = float(bond_table[key2])
            hh = float(bond_table[key3])
            theta = 2.0 * math.asin(hh/2/oh1)
            theta = math.degrees(theta)
            constraints.append((a.index, hs[0], hs[1], 0, 0, 'HOH', theta, oh1, oh2, 0.0, 0.0))

    s = """
    ffio_constraints[%d] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      i_ffio_al
      i_ffio_am
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      r_ffio_c3
      r_ffio_c4
      r_ffio_c5
      :::\n""" % len(constraints)

    i = 1
    for a in constraints:
        s += '      %d  %d  %d  %d  %d  %d  %s  %f  %f  %f  %f  %f\n' % ( i, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10])
        i += 1
    s += '      :::\n'
    s += '    }\n'
    return s

    
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option( '-p', type='str', dest='prmtop_fname', default='',
                       help = 'prmtop file')
    parser.add_option( '-o', type='str', dest='mae_fname', default='',
                       help = 'mae file')
    parser.add_option( '-c', type='str', dest='prmcrd_fname', default='',
                       help = 'prmcrd file')
    opts, args = parser.parse_args()

    blocks = parsePrmtop(opts.prmtop_fname)
    print("Converting crd")
    s1 = convertCrd2Mae(opts.prmcrd_fname, opts.mae_fname)
    print("Omitting top and constraints")
    
    f = open(opts.mae_fname, 'w')
    f.write(s1[:-4])
    f.write('  }\n')
    f.close()

    
