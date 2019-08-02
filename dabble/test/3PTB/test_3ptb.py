# Tests writing amber format with amber parameters
# Special isopeptide bond between two residues
import os
from vmd import atomsel, molecule

dir = os.path.dirname(__file__)
#==============================================================================

def test_multiple_insertion_codes(tmpdir):
    from dabble.param import CharmmWriter

    # Generate the file
    p = str(tmpdir)
    molid = molecule.load("mae", os.path.join(dir, "3PTB_1lig_prepped.mae"))
    w = CharmmWriter(molid=molid, tmp_dir=p, forcefield="charmm", hmr=False,
                     override_defaults=False)
    w.write(os.path.join(p, "test"))

    # Load the output file and start checking it
    m2 = molecule.load("psf", os.path.join(p, "test.psf"),
                       "pdb", os.path.join(p, "test.pdb"))
    molecule.set_top(m2)

    # Check the entire protein is there and connected
    assert len(atomsel("protein or resname ACE NMA NME")) == 3229
    assert len(set(atomsel("protein or resname ACE NMA NME").fragment)) == 1

    # Check the calcium ion is present
    assert atomsel("element Ca").resname == ["CAL"]

    # Check residues with insertion codes
    assert set(atomsel("resid 184").resname) == set(["GLY", "TYR"])
    assert set(atomsel("resid 188").resname) == set(["GLY", "LYS"])
    assert set(atomsel("resid 221").resname) == set(["ALA", "GLN"])
    assert set(atomsel("resid 245").resname) == set(["ASN", "NMA"])
    assert set(atomsel("all").insertion) == set([" ", "A"])

    # Check the ligand is there
    assert len(atomsel("resname BAMI")) == 18

#==============================================================================
