# Tests covalently bonded ligand
import pytest
import subprocess, os

dir = os.path.dirname(__file__) + "/"

#==============================================================================

#def test_multiligand_building(tmpdir):
#    from Dabble import DabbleBuilder
#
#    p = str(tmpdir.mkdir("covligand_build"))
#    print p+"/test.mae"
#    filename =  dir + "rho_arr_CYP_prepped_aligned_dowsered.mae"
#    b = DabbleBuilder(solute_filename=filename, output_filename=p+"/test.mae",
#                      xy_buf=10., wat_buffer=5., overwrite=True, tmp_dir=p)
#    b.write()
#    #resout, reserr = capfd.readouterr()
#    subprocess.check_call(["diff","-q", dir + "test_covligand_correct.mae", p+"/test.mae"])

#==============================================================================

def test_multiligand_renaming(tmpdir):
    """
    Tests covalently bound ligand parameterization with the graph. One
    hydrogen atom should be renamed in the palmitoylcysteine.
    Also tests for detection of phosphorylations on amino acids.
    """
    from DabbleParam import CharmmWriter
    import vmd, molecule

    p = str(tmpdir.mkdir("multiligand_rename"))
    molid = molecule.load("mae", dir + "rho_arr_CYP_prepped_aligned_dowsered.mae")
    w = CharmmWriter(tmp_dir=p, molid=molid, lipid_sel="lipid",
                     extra_topos=[dir + "CYP_v1.str"])
    w.write(p+"/test")
    subprocess.check_call(["diff", "-q", "--ignore-matching-lines=REMARKS",
                           "--ignore-matching-lines=NTITLE",
                           dir + "test_renamed_correct.psf",
                           p+"/test.psf"])

