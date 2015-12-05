# Tests multiple ligands
import pytest
import subprocess


def test_multiligand_building(capfd, tmpdir):
    from Dabble import DabbleBuilder

    p = str(tmpdir.mkdir("multiligand_build"))
    print p+"/test.mae"
    b = DabbleBuilder(solute_filename="B2AR_10ALPs.mae", output_filename=p+"/test.mae",
                      xy_buf=10., wat_buffer=5., overwrite=True, tmp_dir=p)
    b.write()
    #resout, reserr = capfd.readouterr()
    subprocess.check_call(["diff","-q","test_multiligand_correct.mae", p+"/test.mae"])

def test_multiligand_parameterizing(capfd, tmpdir):
    from DabbleParam import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_parameterize"))
    molid = molecule.load("mae", "test_multiligand_correct.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=["alprenolol.rtf"])
    w.write(p+"/test")
    subprocess.check_call(["diff","-q","test_multiligand_correct.psf", p+"/test.psf"])
    subprocess.check_call(["diff","-q","test_multiligand_correct.pdb", p+"/test.pdb"])

def test_multiligand_chamber(capfd, tmpdir):
    from ParmedTools import chamber, parmout
    from chemistry.amber import AmberParm

    p = str(tmpdir.mkdir("multiligand_chamber"))
    parm = AmberParm()
    action = chamber(parm, "-crd test_multiligand_correct.pdb -psf test_multiligand_correct.psf"
                     " -top alprenolol.rtf -par alprenolol.par")
    action.execute()
    write = parmout(action.parm, "test.prmtop test.inpcrd")
    write.execute()
    subprocess.check_call(["diff","-q","test_multiligand_correct.prmtop", p+"/test.prmtop"])
    subprocess.check_call(["diff","-q","test_multiligand_correct.inpcrd", p+"/test.inpcrd"])




