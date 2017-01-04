# Tests molecule matcher
# Im attempting test driven developmetn

import pytest
import subprocess
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
    from Dabble.param import AmberMatcher

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
    from an off file
    """
    from Dabble.param import AmberMatcher
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
    assert("-" not in ala.nodes())

#==============================================================================

def test_residue_renaming(tmpdir):
    """
    Reads a mol2 of LSD with incorrect names as a VMD graph.
    Then, loads a lib file for LSD with correct names, and gaff
    atom types. Should ideally match them up correctly.
    """
    import vmd, molecule
    from atomsel import atomsel
    from Dabble.param import AmberMatcher, MoleculeMatcher

    # Generate a leaprc file with correct paths
    tmpdir = str(tmpdir)
    filename = os.path.join(tmpdir, "leaprc.lsd")
    fh = open(filename, 'w')
    fh.write("source %s\n" % os.path.join(dir, "leaprc.gaff"))
    fh.write("loadOff %s\n" % os.path.join(dir, "lsd.lib"))
    fh.close()

    molid = molecule.load("mae", os.path.join(dir, "lsd_prot.mae"))
    g = AmberMatcher([filename])
    resnaem, mdict = g.get_names(atomsel())

    # Check matching. Ignore hydrogens since those can vary
    assert set(resnaem.values()) == set(["LIG"])
    nohdict = dict((k,v) for k,v in mdict.iteritems() if v[0] != "H")
    assert nohdict == {0: 'C13', 1: 'N2', 2: 'C12', 3: 'C11', 4: 'C15', 5: 'C7',
                       6: 'C3', 7: 'C8', 8: 'N1', 9: 'C1', 10: 'C4', 11: 'C5',
                       12: 'C2', 13: 'C6', 14: 'C9', 15: 'C10', 16: 'C14',
                       17: 'C16', 18: 'O1', 19: 'N3', 20: 'C19', 21: 'C20',
                       22: 'C17', 23: 'C18'}

#==============================================================================

