# Tests covalently bonded ligand
import pytest
import subprocess, os
try:
    import vmd, molecule
    from atomsel import atomsel
except ImportError:
    from vmd import atomsel, molecule
    atomsel = atomsel.atomsel


dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_water_box(tmpdir):
    """
    Tests building in a water box only
    """
    from Dabble import DabbleBuilder
    from Dabble import molutils

    # Build a system with well defined dimensions
    p = str(tmpdir.mkdir("po4_hmr"))
    filename =  os.path.join(dir, "rho_test.mae")
    b = DabbleBuilder(solute_filename=filename,
                      output_filename=os.path.join(p, "test.mae"),
                      membrane_system="TIP3",
                      user_x=40.0, user_y=35.0, user_z = 65.0,
                      tmp_dir=p)
    b.write()

    # Load the built system for checking
    m2 = molecule.load("mae", os.path.join(p, "test.mae"))
    molecule.set_top(m2)

    # Check the system dimensions
    assert molutils.get_system_dimensions(m2) == (40.0, 35.0, 65.0)

#==============================================================================

def test_hmr_param(tmpdir):
    """
    Tests phosphorylations on Ser
    Also checks HMR
    """
    from Dabble.param import AmberWriter
    try:
        import vmd, molecule
    except ImportError:
        from vmd import  molecule

    # Build the system with HMR
    p = str(tmpdir.mkdir("hmr_param"))
    molid = molecule.load("mae", os.path.join(dir, "rho_test.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, hmr=True,
                    extra_topos=[], extra_params=[])
    w.write(os.path.join(p, "test"))

    # Load the built system for checking
    # Use psf so we can check HMR was done correctly
    m2 = molecule.load("psf", os.path.join(p, "test.psf"))
    m3 = molecule.load("parm7", os.path.join(p, "test.prmtop"))

    # Check the protein is there with correct capping groups
    assert len(atomsel("protein", m2)) == 299
    assert len(atomsel("protein", m3)) == 299
    assert len(set(atomsel("resname ACE NMA", m2).get("resid"))) == 2
    assert len(set(atomsel("resname ACE NMA", m3).get("resid"))) == 2

    # Check for 3 phosphoserines and no regular serines
    molecule.set_top(m2)
    assert atomsel("resname SER").get("name").count("P") == 3
    assert not atomsel("resname SER and not same residue as name P")

    # Check for one phosphothreonine
    assert len(set(atomsel("resname THR").get("resid"))) == 3
    assert len(atomsel("resname THR and name P")) == 1

    # Check HMR was done correctly
    minmasspre = min(atomsel("all", m2).get("mass"))
    minmasspost = min(atomsel("all", m3).get("mass"))
    assert abs(minmasspre - 1.008) < 0.0001
    assert abs(minmasspost - 1.008) > 0.0001

#==============================================================================

