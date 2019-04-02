# Tests multiple ligands
import os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_multiligand_building(tmpdir):
    """
    Solvates and membranes a system with multiple ligands """
    from vmd import atomsel, molecule
    from dabble import DabbleBuilder

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

def test_multiligand_parameterizing(tmpdir):
    """
    Checks the parameterization of a system with multiple ligands
    """
    from vmd import atomsel, molecule
    from dabble.param import CharmmWriter

    # Parameterize the multiple ligand system with charmm parameters
    p = str(tmpdir.mkdir("multiligand_parameterize"))
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[os.path.join(dir, "alprenolol.rtf")])
    w.write(os.path.join(p, "test"))

    # Load the created file
    m2 = molecule.load("psf", os.path.join(p, "test.psf"),
                       "pdb", os.path.join(p, "test.pdb"))
    molecule.set_top(m2)

    # Check the system is intact
    assert len(set(atomsel("protein").resid)) == 282
    assert len(set(atomsel("resname ACE NMA").resid)) == 4
    assert len(atomsel("water")) == 32106
    assert len(atomsel("lipid")) == 12194

    # Check for the correct number of alprenolols
    assert len(atomsel("resname ALP")) == 420
    assert set(atomsel("resname ALP").resid) == set(range(1, 11))

    # Check coordinates are unique (not zero)
    # It's 419 because two of them have the same x
    assert len(atomsel("resname ALP")) == 420
    assert len(set(atomsel("resname ALP").x)) == 419

    molecule.delete(m2)

#==============================================================================

def test_multiligand_amber(tmpdir):
    """
    Checks that multiple ligands work with AMBER.
    Should rename some atoms in matching
    """
    from vmd import atomsel, molecule
    from dabble.param import AmberWriter

    # Parameterize the system. One atom name will be changed.
    p = str(tmpdir.mkdir("multiligand_amber"))
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, forcefield="amber",
                    extra_topos=[os.path.join(dir, "alp.off")],
                    extra_params=[os.path.join(dir, "alp.frcmod")])
    w.write(os.path.join(p, "test"))

    # Load the built system and see if it works
    m2 = molecule.load("parm7", os.path.join(p, "test.prmtop"),
                       "rst7", os.path.join(p, "test.inpcrd"))
    molecule.set_top(m2)

    # Check results
    assert len(set(atomsel("protein").resid)) == 282
    assert len(set(atomsel("resname ACE NME").resid)) == 4
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

def test_multiligand_chamber(tmpdir):
    """
    Checks that multiple ligands work with AMBER.
    Should rename some atoms in matching
    """
    from vmd import atomsel, molecule
    from dabble.param import AmberWriter

    # Parameterize the system. One atom name will be changed.
    p = str(tmpdir.mkdir("multiligand_amber"))
    molid = molecule.load("mae", os.path.join(dir, "test_multiligand_correct.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, forcefield="charmm",
                    extra_topos=[os.path.join(dir, "alprenolol.rtf")],
                    extra_params=[os.path.join(dir, "alprenolol.prm")])
    w.write(os.path.join(p, "test"))

    # Load result
    m2 = molecule.load("psf", os.path.join(p, "test.psf"),
                       "pdb", os.path.join(p, "test.pdb"))
    molecule.set_top(m2)

    # Check results
    assert len(set(atomsel("protein").resid)) == 282
    assert len(set(atomsel("resname ACE NMA").resid)) == 4
    assert len(atomsel("water")) == 32106
    assert len(atomsel("same fragment as lipid")) == 12194

    # Check for the corrrect number of alprenolols
    assert len(atomsel("resname ALP")) == 420
    assert len(set(atomsel("resname ALP").resid)) == 10
    assert "O1" in set(atomsel("resname ALP").name)
    assert "OX" not in set(atomsel("resname ALP").name)

    # Check coordinates are unique (not zero)
    # Two atoms have the same X but that is okay
    assert len(set(atomsel("resname ALP").x)) == 419

    molecule.delete(m2)

#==============================================================================

