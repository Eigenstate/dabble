# Tests writing amber format with amber parameters
# Special isopeptide bond between two residues
import os
import pytest
from vmd import atomsel, molecule

dir = os.path.dirname(__file__)
#==============================================================================

def test_amber_custom_residues(tmpdir):
    from dabble.param import AmberWriter

    # Generate the file
    p = str(tmpdir.mkdir("3nob_custom"))
    molid = molecule.load("mae", os.path.join(dir, "prepped.mae"))
    w = AmberWriter(molid, tmp_dir=p, forcefield="amber", hmr=False,
                    extra_topos=[os.path.join(dir, "glx.off"),
                                 os.path.join(dir, "lyx.off")],
                    extra_params=[os.path.join(dir, "join.frcmod"),
                                  os.path.join(dir, "analogies.frcmod")],
                    override_defaults=False)
    w.write(os.path.join(p, "test"))

    # Load the output file and start checking it
    m2 = molecule.load("parm7", os.path.join(p, "test.prmtop"),
                       "rst7", os.path.join(p, "test.inpcrd"))
    molecule.set_top(m2)

    # Check the two custom residues are present
    assert len(atomsel("resname GLX")) == 7
    assert len(atomsel("resname LYX")) == 20

    # Check the custom residues have gaff2 atom types
    assert "n" in atomsel("resname LYX").type
    assert "n2" in atomsel("resname GLX").type

    # Check the normal residues have ff14SB atom types
    assert "N" in atomsel("resname LYS").type
    assert "N" in atomsel("resname GLY").type

    # Check that the isopeptide bond is there
    lybonds = []
    for x in atomsel("resname LYX").bonds:
        lybonds.extend(x)
    assert any(x in lybonds for x in atomsel("resname GLX").index)

#==============================================================================

def test_gromacs_amber(tmpdir):
    from dabble.param import GromacsWriter

    p = str(tmpdir.mkdir("3nob_gromacs"))
    molid = molecule.load("mae", os.path.join(dir, "prepped.mae"))
    w = GromacsWriter(molid, tmp_dir=p,
                      forcefield="amber",
                      extra_topos=[os.path.join(dir, "glx.off"),
                                   os.path.join(dir, "lyx.off")],
                      extra_params=[os.path.join(dir, "join.frcmod"),
                                    os.path.join(dir, "analogies.frcmod")],
                      override_defaults=False)
    w.write(os.path.join(p, "test"))

    # Load and check the output file
    m2 = molecule.load("gro", os.path.join(p, "test.gro"))
    molecule.set_top(m2)

    # Check the two custom residues are present
    assert len(atomsel("resname GLX")) == 7
    assert len(atomsel("resname LYX")) == 20

    # Check that the isopeptide bond is there
    lybonds = []
    for x in atomsel("resname LYX").bonds:
        lybonds.extend(x)
    assert any(x in lybonds for x in atomsel("resname GLX").index)

    # Check the parameterized file with parmed API
    from parmed.formats.registry import load_file
    f = load_file(os.path.join(p, "test.top"))

    assert f.residues[153].name == "GLX"
    assert "n2" in [_.type for _ in f.residues[153].atoms]
    assert "n2" not in [_.type for _ in f.residues[152].atoms]
    assert "N" in [_.type for _ in f.residues[152].atoms]

#==============================================================================
