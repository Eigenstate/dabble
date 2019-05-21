# Tests multiple ligands, all forcefield combinations
import os
from vmd import atomsel, molecule
import pytest

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def check_result(format, outdir):
    if format == "prmtop":
        m2 = molecule.load("parm7", os.path.join(outdir, "test.prmtop"),
                           "rst7", os.path.join(outdir, "test.inpcrd"))
    elif format == "psf":
        m2 = molecule.load("psf", os.path.join(outdir, "test.psf"),
                           "pdb", os.path.join(outdir, "test.pdb"))
    elif format == "gro":
        m2 = molecule.load("gro", os.path.join(outdir, "test.gro"))

    molecule.set_top(m2)

    assert len(set(atomsel("protein").resid)) == 282
    assert len(set(atomsel("resname ACE NMA NME").resid)) == 4
    assert len(atomsel("water")) == 32106
    assert len(atomsel("same fragment as lipid")) == 12194

    # Check for the corrrect number of alprenolols
    assert len(atomsel("resname ALP")) == 420
    assert len(set(atomsel("resname ALP").resid)) == 10
    assert "OX" in set(atomsel("resname ALP").name)
    assert "O1" not in set(atomsel("resname ALP").name)

    # Check coordinates are unique (not zero)
    # It's 419 because two of them have the same x
    assert len(atomsel("resname ALP")) == 420
    assert len(set(atomsel("resname ALP").x)) == 419

    molecule.delete(m2)

#==============================================================================

@pytest.mark.skip(reason="debug")
def test_multiligand_building(tmpdir):
    """
    Solvates and membranes a system with multiple ligands """
    from dabble import DabbleBuilder
    p = str(tmpdir)

    # Build the system
    b = DabbleBuilder(solute_filename=os.path.join(dir, "B2AR_10ALPs.mae"),
                      output_filename=os.path.join(p, "test.mae"),
                      xy_buf=10., wat_buffer=10., overwrite=True, tmp_dir=p)
    b.write()

    # Load the built system
    m2 = molecule.load("mae", os.path.join(p, "test.mae"))
    molecule.set_top(m2)

    # Check the system dimensions are correct
    solsel = atomsel("protein or resname ACE NMA ALP")
    prot_x = max(solsel.x) - min(solsel.x)
    prot_y = max(solsel.y) - min(solsel.y)
    prot_z = max(solsel.z) - min(solsel.z)
    assert abs(molecule.get_periodic(m2)["a"] - prot_x - 20.0) < 0.1
    assert abs(molecule.get_periodic(m2)["b"] - prot_y - 20.0) < 0.1
    assert abs(molecule.get_periodic(m2)["c"] - prot_z - 20.0) < 0.1

    # Check all the alprenolols are there
    assert len(set(atomsel("resname ALP").resid)) == 10
    assert len(atomsel("resname ALP")) == 420

    molecule.delete(m2)

#==============================================================================

@pytest.mark.skip(reason="debug")
def test_multiligand_psf_charmm(tmpdir):
    """
    Checks the parameterization of a system with multiple ligands
    """
    from dabble.param import CharmmWriter
    p = str(tmpdir)

    # Parameterize the multiple ligand system with charmm parameters
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[os.path.join(dir, "alprenolol.rtf")])
    w.write(os.path.join(p, "test"))

    # Check result
    check_result("psf", p)

#==============================================================================

@pytest.mark.skip(reason="debug")
def test_multiligand_prmtop_amber(tmpdir):
    """
    Checks that multiple ligands work with AMBER.
    Should rename some atoms in matching
    """
    from dabble.param import AmberWriter
    p = str(tmpdir)

    # Parameterize the system. One atom name will be changed.
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, forcefield="amber",
                    extra_topos=[os.path.join(dir, "alp.off")],
                    extra_params=[os.path.join(dir, "alp.frcmod")])
    w.write(os.path.join(p, "test"))

    # Load the built system and see if it works
    check_result("prmtop", p)

#==============================================================================

@pytest.mark.skip(reason="debug")
def test_multiligand_prmtop_charmm(tmpdir):
    """
    Checks that multiple ligands work with AMBER.
    Should rename some atoms in matching
    """
    from dabble.param import AmberWriter
    p = str(tmpdir)

    # Parameterize the system. One atom name will be changed.
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, forcefield="charmm",
                    extra_topos=[os.path.join(dir, "alprenolol.rtf")],
                    extra_params=[os.path.join(dir, "alprenolol.prm")])
    w.write(os.path.join(p, "test"))

    # Check results
    check_result("psf", p)

#==============================================================================

@pytest.mark.skip(reason="broken")
def test_multiligand_psf_amber(tmpdir):
    """
    Checks multiple ligands, AMBER parameters, CHARMM format
    """
    # TODO: lipids don't work with this. How am I matching them?
    # TODO: lipid names aren't actually being checked...
    from dabble.param import CharmmWriter
    p = str(tmpdir)

    # AMBER params not supported for lipids with CHARMM output
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid, forcefield="amber",
                     extra_topos=[os.path.join(dir, "alp.off")],
                     extra_params=[os.path.join(dir, "alp.frcmod")])
    with pytest.raises(ValueError):
        w.write(os.path.join(p, "test"))

    # This fails because linkages not being parsed?
    molecule.delete(molid)
    molid = molecule.load("mae", os.path.join(dir, "B2AR_10ALPs.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid, forcefield="amber",
                     extra_topos=[os.path.join(dir, "alp.off")],
                     extra_params=[os.path.join(dir, "alp.frcmod")])
    w.write(os.path.join(p, "test"))

    check_result("prmtop", p)

#==============================================================================

@pytest.mark.skip(reason="augh")
def test_multiligand_gro_amber(tmpdir):
    """
    AMBER parameters, GROMACS format
    """
    from dabble.param import GromacsWriter
    p = str(tmpdir)

    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = GromacsWriter(tmp_dir=p, molid=molid, forcefield="amber",
                      extra_topos=[os.path.join(dir, "alp.off")],
                      extra_params=[os.path.join(dir, "alp.frcmod")])
    w.write(os.path.join(p, "test"))

    check_result("gro", p)

#==============================================================================

def test_multiligand_gro_charmm(tmpdir):
    """
    CHARMM parameters, GROMACS format
    """
    from dabble.param import GromacsWriter
    p = str(tmpdir)

    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = GromacsWriter(tmp_dir=p, molid=molid, forcefield="charmm",
                      extra_topos=[os.path.join(dir, "alprenolol.rtf")],
                      extra_params=[os.path.join(dir, "alprenolol.prm")])
    # TODO: this is bonding the lipid to the protein
    # probably some kind of distance based constraint?
    w.write(os.path.join(p, "test"))

    check_result("gro", p)

#==============================================================================
