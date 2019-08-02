"""
Tests several water models
"""
import os
import pytest
from vmd import molecule, atomsel

dir = os.path.dirname(__file__) + "/"

#==============================================================================

def general_water(tmpdir, writerclass, forcefield, water_model):
    """
    Tests water generally
    """
    m = molecule.load("pdb", os.path.join(dir, "tip3.pdb"))
    p = str(tmpdir)

    w = writerclass(molid=m, tmp_dir=p,
                    forcefield=forcefield,
                    water_model=water_model)
    w.write(os.path.join(p, "test"))

    if os.path.isfile(os.path.join(p, "test.prmtop")):
        tm = molecule.load("parm7", os.path.join(p, "test.prmtop"),
                           "rst7", os.path.join(p, "test.inpcrd"))
    else:
        tm = molecule.load("psf", os.path.join(p, "test.psf"),
                           "pdb", os.path.join(p, "test.pdb"))

    molecule.set_top(tm)
    watsel = "water or resname TIP3 TP4E TP4 SPC SPCE"
    assert len(set(atomsel(watsel).residue)) == 1
    assert len(set(atomsel(watsel).resid)) == 1

    if water_model == "tip4e":
        assert len(atomsel()) == 4
    else: # SPC or TIP3
        assert len(atomsel()) == 3

#==============================================================================

def test_tip3_amber(tmpdir):
    from dabble.param import AmberWriter
    general_water(tmpdir, AmberWriter, "amber", "tip3")

#==============================================================================

def test_tip3_charmm(tmpdir):
    from dabble.param import CharmmWriter
    general_water(tmpdir, CharmmWriter, "charmm", "tip3")

#==============================================================================

@pytest.mark.skip(reason="Parmed bug")
def test_tip4e_chamber(tmpdir):
    from dabble.param import AmberWriter
    general_water(tmpdir, AmberWriter, "charmm", "tip4e")

#==============================================================================

def test_tip4e_amber(tmpdir):
    from dabble.param import AmberWriter
    general_water(tmpdir, AmberWriter, "amber", "tip4e")

#==============================================================================

def test_tip4e_charmm(tmpdir):
    from dabble.param import CharmmWriter
    general_water(tmpdir, CharmmWriter, "charmm", "tip4e")

#==============================================================================

def test_spce_charmm(tmpdir):
    from dabble.param import CharmmWriter
    general_water(tmpdir, CharmmWriter, "charmm", "spce")

#==============================================================================

def test_spce_amber(tmpdir):
    from dabble.param import AmberWriter
    general_water(tmpdir, AmberWriter, "amber", "spce")

#==============================================================================

def test_spce_chamber(tmpdir):
    from dabble.param import AmberWriter
    general_water(tmpdir, AmberWriter, "charmm", "spce")

#==============================================================================
