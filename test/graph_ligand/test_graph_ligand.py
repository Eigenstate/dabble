# Tests the MoleculeMatcher functionality by readin in the 
# mae and str file, and matching up atoms

import pytest
import subprocess
import os

import logging
logging.basicConfig()

dir = os.path.dirname(__file__) + "/"

def test_read_str(capfd, tmpdir):
    """
    Tests if a str file can be properly read in
    """
    from DabbleParam import CharmmMatcher
    import networkx as nx

    g = CharmmMatcher([dir+"lsd_prot_trunc.str", dir+"masses.rtf"])
    nx.write_dot(g.known_res["LSD"], str(tmpdir)+"/test.dot")
    subprocess.check_call(["diff", "-q", dir+"correct_str.dot", str(tmpdir)+"/test.dot"])

def test_read_mol(capfd, tmpdir):
    import vmd, molecule
    import networkx as nx
    from atomsel import atomsel
    from DabbleParam import MoleculeMatcher

    molid = molecule.load("mae", dir+"lsd_prot.mae")
    rgraph, dump = MoleculeMatcher.parse_vmd_graph(atomsel())
    nx.write_dot(rgraph, str(tmpdir)+"/test2.dot")
    subprocess.check_call(["diff", "-q", dir+"correct_mae.dot", str(tmpdir)+"/test2.dot"])

def test_compare_mol():
    import vmd, molecule
    from atomsel import atomsel
    from DabbleParam import CharmmMatcher

    molid = molecule.load("mae", dir+"lsd_prot.mae")
    g = CharmmMatcher([dir+"lsd_prot_trunc.str", dir+"masses.rtf"])
    (name, mdict) = g.get_names(atomsel())
    correctnames = {'C19': 20, 'C18': 23, 'H132': 24, 'C13': 0, 'C12': 2, 'C11': 3, 'C10': 15, 'C17': 22, 'C16': 17, 'O1': 18, 'C14': 16, 'H192': 37, 'H202': 39, 'H123': 26, 'H122': 27, 'H121': 28, 'H203': 40, 'C9': 14, 'C8': 7, 'C15': 4, 'H191': 38, 'C3': 6, 'C2': 12, 'C1': 9, 'N20': 1, 'C7': 5, 'C6': 13, 'C5': 11, 'C4': 10, 'C20': 21, 'H10': 35, 'H11': 29, 'H14': 36, 'H171': 42, 'H172': 43, 'N1': 8, 'N3': 19, 'H152': 30, 'H151': 47, 'H8': 31, 'H2': 34, 'H181': 44, 'H201': 41, 'H131': 25, 'H4': 32, 'H5': 33, 'HN1': 48, 'HN2': 49, 'H183': 45, 'H182': 46}
    names = mdict.next()
    assert(name == "LSD")
    assert(names == correctnames)

def test_patches():
    import vmd, molecule
    from atomsel import atomsel
    from DabbleParam import CharmmMatcher 
    from pkg_resources import resource_filename

    molid = molecule.load("mae", dir+"phosphoserine.mae")
    g = CharmmMatcher([dir+"phosphoserine.str"])
    (name, patch, mdict) = g.get_patches(atomsel("resname SEP"))
    correctnames = {'C': 20, '-C': 5, 'O1P': 24, 'CB': 18, 'CA': 17, '+N': 30, 'O': 21, 'N': 16, 'P': 22, 'HN': 26, 'O2P': 25, 'HA': 27, 'OG': 19, 'OT': 23, 'HB1': 28, 'HB2': 29}
    assert(name == "SER")
    assert(patch == "PSEP")
    assert(mdict.next() == correctnames)

def test_protein(tmpdir):
    import vmd, molecule
    import networkx as nx
    from DabbleParam import CharmmMatcher 
    from pkg_resources import resource_filename

    g = CharmmMatcher([resource_filename("DabbleParam", "charmm_parameters/top_all36_prot.rtf")])
    nx.write_dot(g.known_res["TYR"], str(tmpdir)+"/tyr.dot")
    subprocess.check_call(["diff", "-q", dir+"correct_tyr.dot", str(tmpdir)+"/tyr.dot"])
    
