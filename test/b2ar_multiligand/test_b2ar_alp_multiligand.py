# Tests multiple ligands
import pytest
import subprocess

#==============================================================================

def test_multiligand_building(tmpdir):
    from Dabble import DabbleBuilder

    p = str(tmpdir.mkdir("multiligand_build"))
    print p+"/test.mae"
    b = DabbleBuilder(solute_filename="B2AR_10ALPs.mae", output_filename=p+"/test.mae",
                      xy_buf=10., wat_buffer=5., overwrite=True, tmp_dir=p)
    b.write()
    #resout, reserr = capfd.readouterr()
    subprocess.check_call(["diff","-q","test_multiligand_correct.mae", p+"/test.mae"])

#==============================================================================

def test_multiligand_parameterizing(tmpdir):
    from DabbleParam import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_parameterize"))
    molid = molecule.load("mae", "test_multiligand_correct.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=["alprenolol.rtf"])
    w.write(p+"/test")
    subprocess.check_call(["diff","-q","test_multiligand_correct.psf", p+"/test.psf"])
    subprocess.check_call(["diff","-q","test_multiligand_correct.pdb", p+"/test.pdb"])

#==============================================================================

def test_multiligand_chamber(tmpdir):
    from DabbleParam import AmberWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_chamber"))
    molid = molecule.load("mae", "test_multiligand_correct.mae")
    w = AmberWriter(molid=molid, tmp_dir=p, extra_topos=["alprenolol.rtf"],
                    extra_params=["alprenolol.prm"])
    w.write(p+"/test")
    subprocess.check_call(["diff","-q","test_multiligand_correct.prmtop", p+"/test.prmtop"])
    subprocess.check_call(["diff","-q","test_multiligand_correct.inpcrd", p+"/test.inpcrd"])

#==============================================================================

def test_multiligand_renaming(tmpdir):
    from DabbleParam import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae", "B2AR_10ALPs_renamed.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=["alprenolol.rtf"])
    w.write(p+"/test")
    subprocess.check_call(["diff","-q","test_renamed_correct.psf", p+"/test.psf"])
    subprocess.check_call(["diff","-q","test_renamed_correct.pdb", p+"/test.pdb"])


