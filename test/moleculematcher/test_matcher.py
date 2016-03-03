# Tests molecule matcher
# Im attempting test driven developmetn

import pytest
import subprocess
import os

import logging
logging.basicConfig()

dir = os.path.dirname(__file__)

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

def test_atom_naming(tmpdir):
    """
    Tests if atom names are read in correctly from a leaprc file
    """
    from DabbleParam import AmberMatcher

    # Generate a leaprc file with correct paths
    a = AmberMatcher([create_leaprc(str(tmpdir), "ala.lib")])
    assert(len(a.nodenames) == 80)
    assert("CT" in a.nodenames.keys())
    assert("C" in a.nodenames.values())
    assert("LP" not in a.nodenames.keys())
    assert(a.nodenames.get("CA")=="C")

def test_residue_parsing(tmpdir):
    """
    Checks that a residue is assigned the correct names when read
    from an off file
    """
    from DabbleParam import AmberMatcher
    tmpdir = str(tmpdir)

    a = AmberMatcher([create_leaprc(tmpdir, "ala.lib")])
    # Check residue exists
    assert(a.known_res.get("ALA"))
    ala = a.known_res.get("ALA")
    # Check names and types
    assert(ala.node.get("CB").get("type") == "CT")
    assert(ala.node.get("CB").get("resname") == "HALLO")
    # Check connectivity
    assert(("C","CA") in ala.edges())
    assert(("CB","HB1") in ala.edges())
    assert(("HB2","HB3") not in ala.edges())
    # Check extra bonds
    assert(("C","+") in ala.edges())
    assert("-" not in ala.nodes())

def test_residue_renaming(tmpdir):
    """
    Reads a mol2 of LSD with incorrect names as a VMD graph.
    Then, loads a lib file for LSD with correct names, and gaff
    atom types. Should ideally match them up correctly.
    """
    import vmd, molecule
    from atomsel import atomsel
    from DabbleParam import AmberMatcher, MoleculeMatcher

    # Generate a leaprc file with correct paths
    tmpdir = str(tmpdir)
    filename = os.path.join(tmpdir, "leaprc.lsd")
    fh = open(filename, 'w')
    fh.write("source %s\n" % os.path.join(dir, "leaprc.gaff"))
    fh.write("loadOff %s\n" % os.path.join(dir, "lsd.lib"))
    fh.close()

    molid = molecule.load("mae", os.path.join(dir, "lsd_prot.mae"))
    g = AmberMatcher([filename])
    (resnaem, mdict) = g.get_names(atomsel())

    assert(set(resnaem.values()) == set(["LIG"]))
    assert(mdict=={0: 'C13', 1: 'N2', 2: 'C12', 3: 'C11', 4: 'C15', 5: 'C7', 6: 'C3', 7: 'C8', 8: 'N1', 9: 'C1', 10: 'C4', 11: 'C5', 12: 'C2', 13: 'C6', 14: 'C9', 15: 'C10', 16: 'C14', 17: 'C16', 18: 'O1', 19: 'N3', 20: 'C17', 21: 'C18', 22: 'C19', 23: 'C20', 24: 'H131', 25: 'H132', 26: 'H121', 27: 'H122', 28: 'H123', 29: 'H11', 30: 'H151', 31: 'H8', 32: 'H4', 33: 'H5', 34: 'H2', 35: 'H10', 36: 'H14', 37: 'H171', 38: 'H172', 39: 'H181', 40: 'H182', 41: 'H183', 42: 'H191', 43: 'H192', 44: 'H201', 45: 'H202', 46: 'H203', 47: 'H152', 48: 'HN1', 49: 'HN2'})




