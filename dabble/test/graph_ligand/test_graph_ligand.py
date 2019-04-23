# Tests the MoleculeMatcher functionality by readin in the
# mae and str file, and matching up atoms
import os
import logging
import subprocess
logging.basicConfig()

dir = os.path.dirname(__file__) + "/"
#==============================================================================

#def test_read_str(capfd, tmpdir):
#    """
#    Tests if a str file can be properly read in
#    """
#    from dabble.param import CharmmMatcher
#
#    g = CharmmMatcher([dir+"lsd_prot_trunc.str", dir+"masses.rtf"])
#    g.write_dot(g.known_res["LSD"], str(tmpdir)+"/test.dot")
#    subprocess.check_call(["diff", "-q", dir+"correct_str.dot", str(tmpdir)+"/test.dot"])

#==============================================================================

#def test_read_mol(capfd, tmpdir):
#    from vmd import atomsel, molecule
#
#    from dabble.param import MoleculeMatcher
#
#    molid = molecule.load("mae", dir+"lsd_prot.mae")
#    rgraph, dump = MoleculeMatcher.parse_vmd_graph(atomsel())
#    g.write_dot(rgraph, str(tmpdir)+"/test2.dot")
#    subprocess.check_call(["diff", "-q", dir+"correct_mae.dot", str(tmpdir)+"/test2.dot"])

#==============================================================================

def test_compare_mol():
    from vmd import atomsel, molecule
    from dabble.param import CharmmMatcher

    molid = molecule.load("mae", dir+"lsd_prot.mae")
    g = CharmmMatcher([dir+"lsd_prot_trunc.str", dir+"masses.rtf"])
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

def test_patches():
    from vmd import atomsel, molecule
    from dabble.param import CharmmMatcher

    molid = molecule.load("mae", dir+"phosphoserine.mae")
    g = CharmmMatcher([dir+"phosphoserine.str"])
    (name, patch, mdict) = g.get_patches(atomsel("resname SEP", molid=molid))
    assert name == "SER"
    assert patch == "PSEP"

    # Oxygens are interchangeable so don't compare them
    correct = {5: '-C', 16: 'N', 17: 'CA', 18: 'CB', 19: 'OG', 20: 'C', 21:
               'O', 22: 'P', 26: 'HN', 27: 'HA', 28: 'HB1', 29: 'HB2', 30:
               '+N'}
    for i, n in correct.items():
        assert mdict[i] == n

    # Interchangeable oxygens
    assert "O" in mdict[23]
    assert "O" in mdict[24]
    assert "O" in mdict[25]

#==============================================================================

def test_check_resname():
    import pytest
    from dabble import DabbleError
    from dabble.param import CharmmMatcher

    with pytest.raises(DabbleError):
        _ = CharmmMatcher([dir + "lsd_toolong.str"])

#==============================================================================

#def test_protein(tmpdir):
#    from vmd import molecule
#
#    from dabble.param import CharmmMatcher
#    from pkg_resources import resource_filename
#
#    g = CharmmMatcher([resource_filename("Dabble.param", "charmm_parameters/top_all36_prot.rtf")])
#    g.write_dot(g.known_res["TYR"], str(tmpdir)+"/tyr.dot")
#    subprocess.check_call(["diff", "-q", dir+"correct_tyr.dot", str(tmpdir)+"/tyr.dot"])

#==============================================================================

