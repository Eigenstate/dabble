"""
Tests a variety of potentially problematic disulfides. Disulfides between
two same resids on different chains, disulfides on the same chain, and
unrelated disulfides on different chains.

Also tests how chains are handled when Maestro puts a capping group
with resid -1 on them
"""
import pytest
import subprocess, os
try:
    import vmd, molecule
    from atomsel import atomsel
except ImportError:
    from vmd import molecule, atomsel
    atomsel = atomsel.atomsel

dir = os.path.dirname(__file__) + "/"

def check_correctness(molid):
    """ Verifies molecule is sane """
    molecule.set_top(molid)

    # Check the protein is there with the correct capping groups
    assert len(atomsel("protein")) == 804
    assert len(set(atomsel("all").get("fragment"))) == 2
    assert len(set(atomsel("resname ACE NMA NME").get("residue"))) == 4

    # Check for 6 cysteines, 2 with same resid
    assert len(set(atomsel("resname CYS CYX").get("residue"))) == 6

    # Check connectivity between cysteines is correct
    for res in set(atomsel("resname CYS CYX").get("residue")):
        assert len(atomsel("residue %d" % res)) == 10
        assert len(atomsel("residue %d and name SG" % res)) == 1
        idxs = atomsel("residue %d and name SG" % res).bonds[0]
        assert set(atomsel("index %s"
                           % " ".join(str(i) for i in idxs)).get("name")) \
            == set(["CB", "SG"])

#==============================================================================

def test_relaxin_disulfides_charmm(tmpdir):
    """
    Tests these disulfides with charmm
    """
    from Dabble.param import CharmmWriter
    try:
        import vmd, molecule
    except ImportError:
        from vmd import  molecule

    # Build the system
    p = str(tmpdir.mkdir("charmm_disu"))
    molid = molecule.load("mae", os.path.join(dir, "2MV1_dowserwat.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid)
    w.write(os.path.join(p, "test"))

    # Check the built system
    m2 = molecule.load("psf", os.path.join(p, "test.psf"))
    check_correctness(m2)

#==============================================================================

def test_relaxin_disulfides_amber(tmpdir):
    """
    Tests these disulfides with amber
    """
    from Dabble.param import AmberWriter
    try:
        import vmd, molecule
        from atomsel import atomsel
    except ImportError:
        from vmd import molecule, atomsel
        atomsel = atomsel.atomsel

    # Build the system
    p = str(tmpdir.mkdir("amber_disu"))
    molid = molecule.load("mae", os.path.join(dir, "2MV1_dowserwat.mae"))
    w = AmberWriter(molid=molid, tmp_dir=p, forcefield="amber")
    w.write(os.path.join(p, "test"))

    # Load and sanity check the built system
    m2 = molecule.load("parm7", os.path.join(p, "test.prmtop"))
    check_correctness(m2)


