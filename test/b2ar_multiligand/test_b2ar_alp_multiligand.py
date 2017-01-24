# Tests multiple ligands
import pytest
import subprocess, os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_multiligand_building(tmpdir):
    """
    Solvates and membranes a system with multiple ligands
    """
    try:
        import vmd, molecule
        from atomsel import atomsel
    except ImportError:
        from vmd import atomsel, molecule
        atomsel = atomsel.atomsel

    from Dabble import DabbleBuilder

    # Build the system
    p = str(tmpdir.mkdir("multiligand_build"))
    b = DabbleBuilder(solute_filename=os.path.join(dir, "B2AR_10ALPs.mae"),
                      output_filename=os.path.join(p, "test.mae"),
                      xy_buf=10., wat_buffer=10., overwrite=True, tmp_dir=p)
    b.write()

    # Load the built system
    m2 = molecule.load("mae", os.path.join(p, "test.mae"))
    molecule.set_top(m2)

    # Check the system dimensions are correct
    solsel = atomsel("protein or resname ACE NMA ALP")
    prot_x = max(solsel.get("x")) - min(solsel.get("x"))
    prot_y = max(solsel.get("y")) - min(solsel.get("y"))
    prot_z = max(solsel.get("z")) - min(solsel.get("z"))
    assert(abs(molecule.get_periodic(m2)["a"] - prot_x - 20.0) < 0.1)
    assert(abs(molecule.get_periodic(m2)["b"] - prot_y - 20.0) < 0.1)
    assert(abs(molecule.get_periodic(m2)["c"] - prot_z - 20.0) < 0.1)

    # Check all the alprenolols are there
    assert(len(set(atomsel("resname ALP").get("resid"))) == 10)

#==============================================================================

def test_multiligand_parameterizing(tmpdir):
    """
    Checks the parameterization of a system with multiple ligands
    """
    try:
        import vmd, molecule
        from atomsel import atomsel
    except ImportError:
        from vmd import atomsel, molecule
        atomsel = atomsel.atomsel

    from Dabble.param import CharmmWriter

    # Parameterize the multiple ligand system with charmm parameters
    p = str(tmpdir.mkdir("multiligand_parameterize"))
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[os.path.join(dir, "alprenolol.rtf")])
    w.write(os.path.join(p, "test"))

    # Load the created file
    m2 = molecule.load("psf", os.path.join(p, "test.psf"))
    molecule.set_top(m2)

    # Check the system is intact
    assert(len(set(atomsel("protein").get("resid"))) == 282)
    assert(len(set(atomsel("resname ACE NMA").get("resid"))) == 4)
    assert(len(atomsel("water")) == 32106)
    assert(len(atomsel("lipid")) == 12194)

    # Check for the correct number of alprenolols
    assert(len(atomsel("resname ALP")) == 420)
    assert(set(atomsel("resname ALP").get("resid")) == set(range(1,11)))

#==============================================================================
# Commented out because really slow and chamber has its own tests
#def test_multiligand_chamber(tmpdir):
#    from Dabble.param import AmberWriter
#    import vmd, molecule
#
#    p = str(tmpdir.mkdir("multiligand_chamber"))
#    molid = molecule.load("mae", dir + "test_multiligand_correct.mae")
#    w = AmberWriter(molid=molid, tmp_dir=p, extra_topos=[dir + "alprenolol.rtf"],
#                    extra_params=[dir + "alprenolol.prm"])
#    w.write(p+"/test")
#    subprocess.check_call(["diff", "-q", dir + "test_multiligand_correct.prmtop",
#                           p+"/test.prmtop"])
#    subprocess.check_call(["diff", "-q", dir + "test_multiligand_correct.inpcrd",
#                           p+"/test.inpcrd"])

#==============================================================================

#def test_multiligand_renaming(tmpdir):
#    from Dabble.param import CharmmWriter
#    import vmd, molecule
#
#    p = str(tmpdir.mkdir("multiligand_rename"))
#    molid = molecule.load("mae", dir + "B2AR_10ALPs_renamed.mae")
#    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
#                     extra_topos=[dir + "alprenolol.rtf"])
#    w.write(p+"/test")
#    subprocess.check_call(["diff", "-q", "--ignore-matching-lines=REMARKS",
#                           dir + "test_renamed_correct.psf",
#                           "--ignore-matching-lines=NTITLE",
#                           p+"/test.psf"])
#
#==============================================================================

#def test_amber_writing(tmpdir):
#    """
#    Tests saving a system with no ligand to a file
#    """
#    from Dabble.param import AmberWriter
#    import vmd, molecule
#
#    p = str(tmpdir.mkdir("multiligand_rename"))
#    molid = molecule.load("mae", dir + "B2AR_noalps.mae")
#    w = AmberWriter(molid, tmp_dir=p, forcefield="amber")
#    w.write(os.path.join(p,"test"))

#==============================================================================

#def test_multiligand_chamber(tmpdir):
#    from Dabble.param import AmberWriter
#    import vmd, molecule
#
##    p = str(tmpdir.mkdir("multiligand_rename"))
#    p = tmpdir
#    dir = "./"
#    molid = molecule.load("mae", dir + "test_multiligand_noalps.mae")
#    w = AmberWriter(molid, tmp_dir=p, forcefield="charmm", hmr=True)
#    w.write("test_chamber")
#
#==============================================================================

