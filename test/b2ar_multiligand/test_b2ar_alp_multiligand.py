# Tests multiple ligands
import pytest
import subprocess, os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def test_multiligand_building(tmpdir):
    from Dabble import DabbleBuilder

    p = str(tmpdir.mkdir("multiligand_build"))
    print p+"/test.mae"
    filename =  dir + "B2AR_10ALPs.mae"
    b = DabbleBuilder(solute_filename=filename, output_filename=p+"/test.mae",
                      xy_buf=5., wat_buffer=5., overwrite=True, tmp_dir=p)
    b.write()
    #resout, reserr = capfd.readouterr()
    subprocess.check_call(["diff","-q", dir + "test_multiligand_correct.mae", p+"/test.mae"])

#==============================================================================

def test_multiligand_parameterizing(tmpdir):
    from Dabble.param import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_parameterize"))
    molid = molecule.load("mae", dir + "test_multiligand_correct.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[dir+"alprenolol.rtf"])
    w.write(p+"/test")
    subprocess.check_call(["diff", "-q", "--ignore-matching-lines=REMARKS",
                           "--ignore-matching-lines=NTITLE",
                           dir + "test_multiligand_correct.psf", p+"/test.psf"])

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

def test_multiligand_renaming(tmpdir):
    from Dabble.param import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae", dir + "B2AR_10ALPs_renamed.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[dir + "alprenolol.rtf"])
    w.write(p+"/test")
    subprocess.check_call(["diff", "-q", "--ignore-matching-lines=REMARKS",
                           dir + "test_renamed_correct.psf",
                           "--ignore-matching-lines=NTITLE",
                           p+"/test.psf"])


