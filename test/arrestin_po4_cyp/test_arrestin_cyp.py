# Tests covalently bonded ligand
import os

dir = os.path.dirname(__file__)

#==============================================================================

def test_covalent_ligand_patches(tmpdir):
    """
    Tests covalently bound ligand parameterization with the graph. One
    hydrogen atom should be renamed in the palmitoylcysteine.
    Also tests for detection of phosphorylations on amino acids.
    """
    from vmd import atomsel, molecule
    from Dabble.param import CharmmWriter

    # Parameterize with charmm parameters
    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae",
                          os.path.join(dir,
                                       "rho_arr_CYP_prepped_aligned_dowsered.mae"))

    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[os.path.join(dir, "CYP_v1.str")])
    w.write(os.path.join(p, "test"))

    # Load the result
    m2 = molecule.load("psf", os.path.join(p, "test.psf"))
    molecule.set_top(m2)

    # Sanity check the system was built completely
    # Some resids match because of insertion codes
    assert(len(set(atomsel("protein or resname ACE NMA").get("resid"))) == 697)
    assert(len(set(atomsel("protein or resname ACE NMA").get("residue"))) == 699)
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

    # Check for protonation on Glu134A
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

def test_covalent_ligand_amber(tmpdir):
    """
    Tests covalently bound ligand parameterization.
    Also tests for detection of phosphorylations on amino acids.

    Sanity checks differ here from the charmm once since numbering in
    the prmtop file is different.
    """
    from vmd import atomsel, molecule
    from Dabble.param import AmberWriter

    # Parameterize with charmm parameters
    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae",
                          os.path.join(dir,
                                       "rho_arr_CYP_prepped_aligned_dowsered.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                    forcefield="amber",
                    extra_topos=[os.path.join(dir, "cyp.off"),
                                 os.path.join(os.environ.get("AMBERHOME"),
                                              "dat", "leap", "cmd",
                                              "leaprc.phosaa10")],
                    extra_params=[os.path.join(dir, "cyp.frcmod"),
                                  os.path.join(dir, "analogies.frcmod")]
                   )
    w.write(os.path.join(p, "test"))

    # Load the result
    m2 = molecule.load("parm7", os.path.join(p, "test.prmtop"))
    molecule.set_top(m2)

    # Sanity check the system was built completely
    assert(len(set(atomsel("protein or resname ACE NME").get("resid"))) == 699)
    assert(len(set(atomsel("protein or resname ACE NME").get("residue"))) == 699)
    assert(len(set(atomsel("resname ACE NME").get("resid"))) == 4)
    assert(len(atomsel("water")) == 654)
    assert(len(set(atomsel("all").get("fragment"))) == 220)

    # Check for palmitoylation
    assert(len(atomsel("resname CYP")) == 116)
    assert(set(atomsel("resname CYP").get("resid")) == set([323, 324]))

    # Check for protonation on Asp83
    # Check for no protonation on Asp334
    assert(len(atomsel("resid 84 and resname ASH")) == 13)
    assert(len(atomsel("resid 332 and resname ASP")) == 12)
    assert("OH" in atomsel("resid 84 and resname ASH").get("type"))
    assert("OH" not in atomsel("resid 332 and resname ASP").get("type"))

    # Check for protonation on Glu134A
    assert(len(atomsel("resid 135 and resname GLH")) == 16)
    assert(len(atomsel("resid 248 and resname GLU")) == 15)
    assert("OH" in atomsel("resid 135 and resname GLH").get("type"))
    assert("OH" not in atomsel("resid 332 and resname GLU").get("type"))

    # Check that Ser98 is normal
    assert(len(atomsel("resid 99 and resname SER")) == 11)
    assert("P" not in atomsel("resid 99 and resname SER").get("type"))
    assert("OH" in atomsel("resid 99 and resname SER").get("type"))

    # Check for phosphorylation on Ser334
    assert(len(atomsel("resid 335 and resname SEP")) == 14)
    assert("P" in atomsel("resid 335 and resname SEP").get("type"))
    assert("OH" not in atomsel("resid 335 and resname SEP").get("type"))

    # Check for phosphorylation on Ser338
    assert(len(atomsel("resid 339 and resname SEP")) == 14)
    assert("P" in atomsel("resid 339 and resname SEP").get("type"))
    assert("OH" not in atomsel("resid 339 and resname SEP").get("type"))

    # Check for disulfide bonds
    assert all(len(x)==2 for x in atomsel("resname CYX and element S").bonds)

#==============================================================================
