# Tests covalently bonded ligand
import pytest
import subprocess, os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_water_box(tmpdir):
    """
    Tests building in a water box only
    """
    from Dabble import DabbleBuilder

    p = str(tmpdir.mkdir("po4_hmr"))
#    print p+"/test.mae"
    filename =  dir + "rho_test.mae"
    b = DabbleBuilder(solute_filename=filename, output_filename=p+"/test.mae",
                      membrane_system="TIP3", wat_buffer=5.,
                      overwrite=True, tmp_dir=p)
    b.write()
    #resout, reserr = capfd.readouterr()
    subprocess.check_call(["diff","-q", dir + "test_rho_correct.mae", p+"/test.mae"])

#==============================================================================

def test_hmr_param(tmpdir):
    """
    Tests phosphorylations on Ser, Thr
    Also checks HMR
    """
    from DabbleParam import AmberWriter 
    import vmd, molecule

    p = str(tmpdir.mkdir("hmr_param"))
    molid = molecule.load("mae", dir + "test_rho_correct.mae")
    w = AmberWriter(tmp_dir=p, molid=molid, hmr=True,
                    extra_topos=[], extra_params=[])
    w.write(p+"/test")
    subprocess.check_call(["diff", "-q", "--ignore-matching-lines=VERSION",
                           dir + "test_rho_correct.prmtop",
                           p+"/test.prmtop"])

#==============================================================================

