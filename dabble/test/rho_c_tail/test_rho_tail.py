# Tests covalently bonded ligand
import os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def check_built_system(p):

    from vmd import atomsel, molecule

    # Load the built system for checking
    assert os.path.isfile(os.path.join(p, "test.prmtop"))

    # Use psf so we can check HMR was done correctly
    m3 = molecule.load("parm7", os.path.join(p, "test.prmtop"))

    # Check the protein is there with correct capping groups
    assert len(atomsel("protein", m3)) == 299
    assert len(set(atomsel("resname ACE NMA NME", m3).resid)) == 2

    # Check for 3 phosphoserines and no regular serines
    molecule.set_top(m3)
    assert atomsel("resname SER SEP").name.count("P") == 3
    assert not atomsel("resname SER SEP and not same residue as name P")

    # Check for one phosphothreonine
    assert len(set(atomsel("resname THR TPO").resid)) == 3
    assert len(atomsel("resname THR TPO and name P")) == 1

    # Check HMR was done correctly
    minmasspost = min(atomsel("all", m3).mass)
    assert abs(minmasspost - 1.008) > 0.0001

#==============================================================================

def test_water_box(tmpdir):
    """
    Tests building in a water box only
    """
    from dabble import DabbleBuilder
    from dabble import molutils
    from vmd import molecule

    # Build a system with well defined dimensions
    p = str(tmpdir)
    filename = os.path.join(dir, "rho_test.mae")
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

def test_hmr_param_amber(tmpdir):
    """
    Tests parameterizing with Amber, using phosphoserines
    """
    from dabble.param import AmberWriter
    from vmd import molecule

    # Build system
    p = str(tmpdir)
    molid = molecule.load("mae", os.path.join(dir, "rho_test.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, hmr=True, forcefield="amber",
                    extra_topos=[os.path.join(os.environ.get("AMBERHOME"),
                                              "dat", "leap", "cmd",
                                              "leaprc.phosaa10")])
    w.write(os.path.join(p, "test"))
    check_built_system(p)

#==============================================================================

def test_hmr_param(tmpdir):
    """
    Tests phosphorylations on Ser
    Also checks HMR
    """
    from dabble.param import AmberWriter
    from vmd import molecule

    # Build the system with HMR
    p = str(tmpdir)
    molid = molecule.load("mae", os.path.join(dir, "rho_test.mae"))
    w = AmberWriter(tmp_dir=p, molid=molid, hmr=True, forcefield="charmm")
    w.write(os.path.join(p, "test"))
    check_built_system(p)

#==============================================================================

