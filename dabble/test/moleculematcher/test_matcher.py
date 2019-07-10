# Tests molecule matcher
# Im attempting test driven developmetn

import os

import logging
logging.basicConfig()

dir = os.path.dirname(__file__)

#==============================================================================

def create_leaprc(tmpdir, libfile):
    """
    Creates a leaprc pointing to libfile with atom types
    """
    import shutil
    filename = os.path.abspath(os.path.join(tmpdir, "leaprc.test"))
    shutil.copy2(os.path.join(dir,"leaprc.skelly"), filename)
    a=open(filename, 'a')
    a.write("loadOff %s" % os.path.join(dir,libfile))
    a.close()
    return os.path.abspath(filename)

#==============================================================================

def test_atom_naming(tmpdir):
    """
    Tests if atom names are read in correctly from a leaprc file
    """
    from dabble.param import AmberMatcher

    # Generate a leaprc file with correct paths
    a = AmberMatcher([create_leaprc(str(tmpdir), "ala.lib")])
    assert len(a.nodenames) == 80
    assert "CT" in a.nodenames.keys()
    assert "C" in a.nodenames.values()
    assert "LP" not in a.nodenames.keys()
    assert a.nodenames.get("CA")=="C"

#==============================================================================

def test_residue_parsing(tmpdir):
    """
    Checks that a residue is assigned the correct names when read
    from an off file in AMBER format
    """
    from dabble.param import AmberMatcher
    tmpdir = str(tmpdir)

    a = AmberMatcher([create_leaprc(tmpdir, "ala.lib")])

    # Check residue exists
    assert a.known_res.get("ALA")
    ala = a.known_res.get("ALA")

    # Check names and types
    print(ala.node.keys())
    assert "CB" in [ala.node[x].get("atomname") for x in ala.nodes()]
    cbs = [x for x in ala.nodes() if ala.node[x].get("atomname") == "CB"]
    assert len(cbs) == 1
    assert cbs[0] == "5"
    assert ala.node[cbs[0]].get("type") and ala.node[cbs[0]].get("type") == "CT"
    assert set([ala.node[x].get("resname") for x in ala.nodes()]) == set(["HALLO", None])

    # Check connectivity
    C = [x for x in ala.nodes() if ala.node[x].get("atomname") == "C"][0]
    CA = [x for x in ala.nodes() if ala.node[x].get("atomname") == "CA"][0]
    CB = [x for x in ala.nodes() if ala.node[x].get("atomname") == "CB"][0]
    HB1 = [x for x in ala.nodes() if ala.node[x].get("atomname") == "HB1"][0]
    HB2 = [x for x in ala.nodes() if ala.node[x].get("atomname") == "HB2"][0]
    HB3 = [x for x in ala.nodes() if ala.node[x].get("atomname") == "HB3"][0]

    assert (C, CA) in ala.edges() or (CA, C) in ala.edges()
    assert (CB, HB1) in ala.edges() or (HB1, CB) in ala.edges()
    assert (HB2, HB3) not in ala.edges() and (HB3, HB2) not in ala.edges()

    # Check extra bonds
    assert (C, "+") in ala.edges() or ("+", C) in ala.edges()
    assert "-" not in ala.nodes()

#==============================================================================

def test_residue_renaming(tmpdir):
    """
    Reads a mol2 of LSD with incorrect names as a VMD graph.
    Then, loads a lib file for LSD with correct names, and gaff
    atom types. Tests AMBER matching.

    The LSD definition also has a pseudoatom that needs to be correctly
    ignored in order for matching to occur.
    """
    from vmd import atomsel, molecule
    from dabble.param import AmberMatcher

    # Generate a leaprc file with correct paths
    tmpdir = str(tmpdir)
    filename = os.path.join(tmpdir, "leaprc.lsd")
    fh = open(filename, 'w')
    fh.write("source %s\n" % os.path.join(dir, "leaprc.gaff"))
    fh.write("loadOff %s\n" % os.path.join(dir, "lsd.lib"))
    fh.close()

    molecule.load("mae", os.path.join(dir, "lsd_prot.mae"))
    g = AmberMatcher([filename])
    resnaem, mdict = g.get_names(atomsel())

    assert resnaem is not None

    # Check matching. Ignore hydrogens since those can vary
    assert set(resnaem.values()) == set(["LIG"])
    nohdict = dict((k,v) for k,v in mdict.items() if v[0] != "H")
    assert nohdict == {0: 'C13', 1: 'N2', 2: 'C12', 3: 'C11', 4: 'C15', 5: 'C7',
                       6: 'C3', 7: 'C8', 8: 'N1', 9: 'C1', 10: 'C4', 11: 'C5',
                       12: 'C2', 13: 'C6', 14: 'C9', 15: 'C10', 16: 'C14',
                       17: 'C16', 18: 'O1', 19: 'N3', 20: 'C19', 21: 'C20',
                       22: 'C17', 23: 'C18'}

#==============================================================================

def test_patches():
    """
    Tests molecule matching in CHARMM with patches
    """
    from vmd import atomsel, molecule
    from dabble.param import CharmmMatcher

    molid = molecule.load("mae", os.path.join(dir, "phosphoserine.mae"))
    g = CharmmMatcher([os.path.join(dir, "phosphoserine.str")])
    (name, patch, mdict) = g.get_patches(atomsel("resname SEP", molid=molid))
    assert name == "SER"
    assert patch == "PSEP"

    # Oxygens are interchangeable so don't compare them
    correct = {
        16: 'N', 17: 'CA', 18: 'CB', 19: 'OG', 20: 'C', 21: 'O', 22: 'P',
        26: 'HN', 27: 'HA', 28: 'HB1', 29: 'HB2'
               }
    for i, n in correct.items():
        assert mdict[i] == n

    # Interchangeable oxygens
    assert "O" in mdict[23]
    assert "O" in mdict[24]
    assert "O" in mdict[25]

#==============================================================================

def test_check_resname():
    """
    Tests check for 6 character resname in CHARMM
    """
    import pytest
    from dabble import DabbleError
    from dabble.param import CharmmMatcher

    with pytest.raises(DabbleError):
        _ = CharmmMatcher([os.path.join(dir, "lsd_toolong.str")])

#==============================================================================

def test_compare_mol():
    """
    Tests name matching with CHARMM
    """
    from vmd import atomsel, molecule
    from dabble.param import CharmmMatcher

    molid = molecule.load("mae", os.path.join(dir, "lsd_prot.mae"))
    g = CharmmMatcher([os.path.join(dir, "lsd_prot_trunc.str"),
                       os.path.join(dir, "masses.rtf")])
    (resdict, mdict) = g.get_names(atomsel("all", molid=molid))

    # Don't check diethyls as they're interchangeable
    correct = {0: 'C13', 1: 'N20', 2: 'C12', 3: 'C11', 4: 'C15', 5: 'C7',
               6: 'C3', 7: 'C8', 8: 'N1', 9: 'C1', 10: 'C4', 11: 'C5', 12: 'C2',
               13: 'C6', 14: 'C9', 15: 'C10', 16: 'C14', 17: 'C16', 18: 'O1',
               19: 'N3', 24: 'H131', 25: 'H132', 26: 'H121', 27: 'H122',
               28: 'H123', 29: 'H11', 30: 'H151', 31: 'H8', 32: 'H4',
               33: 'H5', 34: 'H2', 35: 'H10', 36: 'H14', 47: 'H152',
               48: 'HN1', 49: 'HN2'}
    for i, n in correct.items():
        assert mdict[i] == n

    assert "LSD" == resdict

#==============================================================================

def test_gromacs_itp():

    from vmd import atomsel, molecule
    from dabble.param import GromacsMatcher

    molid = molecule.load("mae", os.path.join(dir, "lsd_prot.mae"))
    g = GromacsMatcher(topologies=[os.path.join(dir, "lsd.itp")])

    print(list(g.known_res.keys()))
    assert "GLSD" in g.known_res.keys()
    g.write_dot(g.known_res["GLSD"], "test.dot")
    (resdict, mdict) = g.get_names(atomsel("all", molid))

    correct = {13: 'C12', 6: 'C13', 5: 'C14', 7: 'C15', 8: 'N16', 9: 'C17', 10:
               'C18', 11: 'C19', 12: 'C20', 14: 'C11', 15: 'C10', 16: 'C8', 17:
               'C6', 19: 'N3', 20: 'C2', 21: 'C1', 22: 'C4', 23: 'C5', 18:
               'O7', 36: 'H9', 4: 'C21', 3: 'C22', 29: 'H23', 1: 'N24', 49:
               'H25', 2: 'C26', 0: 'C27', 39: 'H28', 40: 'H29', 41: 'H30', 37:
               'H31', 38: 'H32', 42: 'H33', 43: 'H34', 44: 'H35', 45: 'H36',
               46: 'H37', 35: 'H38', 31: 'H39', 48: 'H40', 32: 'H41', 33:
               'H42', 34: 'H43', 47: 'H44', 30: 'H45', 26: 'H46', 27: 'H47',
               28: 'H48', 24: 'H49', 25: 'H50'}

    assert all(k == "GLSD" for k in resdict.values())
    for i, n in correct.items():
        assert mdict[i] == n

#==============================================================================
