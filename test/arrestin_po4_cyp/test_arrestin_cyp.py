# Tests covalently bonded ligand
import pytest
import subprocess, os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_covalent_ligand_patches(tmpdir):
    """
    Tests covalently bound ligand parameterization with the graph. One
    hydrogen atom should be renamed in the palmitoylcysteine.
    Also tests for detection of phosphorylations on amino acids.
    """
    try:
        import vmd, molecule
        from atomsel import atomsel
    except ImportError:
        from vmd import atomsel, molecule
        atomsel = atomsel.atomsel

    from Dabble.param import CharmmWriter

    # Parameterize with charmm parameters
    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae", dir + "rho_arr_CYP_prepped_aligned_dowsered.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[dir + "CYP_v1.str"])
    w.write(os.path.join(p, "test"))

    # Load the result
    m2 = molecule.load("psf", os.path.join(p, "test.psf"))
    molecule.set_top(m2)

    # Sanity check the system was built completely
    assert(len(set(atomsel("protein or resname ACE NMA").get("resid"))) == 699)
    assert(len(set(atomsel("resname ACE NMA").get("resid"))) == 4)
    assert(len(atomsel("water")) == 654)
    assert(len(set(atomsel("all").get("fragment"))) == 220)

    # Check for palmitoylation
    assert(len(atomsel("resname CYP")) == 116)
    assert(set(atomsel("resname CYP").get("resid")) == set([322, 323]))

    # Check for protonation on Asp83
    assert(len(atomsel("resid 83 and resname ASP")) == 13)
    assert(len(atomsel("resid 331 and resname ASP")) == 12)
    assert("OH1" in atomsel("resid 83 and resname ASP").get("type"))
    assert("OH1" not in atomsel("resid 331 and resname ASP").get("type"))

    # Check for protation on Glu134A
    assert(len(atomsel("resid 134 and resname GLU")) == 16)
    assert(len(atomsel("resid 247 and resname GLU")) == 15)
    assert("OB" in atomsel("resid 134 and resname GLU").get("type"))
    assert("OB" not in atomsel("resid 331 and resname GLU").get("type"))

    # Check that Ser98 is normal
    assert(len(atomsel("resid 98 and resname SER")) == 11)
    assert("P" not in atomsel("resid 98 and resname SER").get("type"))
    assert("OH1" in atomsel("resid 98 and resname SER").get("type"))

    # Check for phosphorylation on Ser334
    assert(len(atomsel("resid 334 and resname SER")) == 14)
    assert("P" in atomsel("resid 334 and resname SER").get("type"))
    assert("OH1" not in atomsel("resid 334 and resname SER").get("type"))

    # Check for phosphorylation on Ser338
    assert(len(atomsel("resid 338 and resname SER")) == 14)
    assert("P" in atomsel("resid 338 and resname SER").get("type"))
    assert("OH1" not in atomsel("resid 338 and resname SER").get("type"))

#==============================================================================
