"""
Tests absolute box size with ignored selection
"""
import os
import pytest

dir = os.path.dirname(__file__) + "/"

#==============================================================================
@pytest.mark.skip(reason="Missing input file")
def test_absolute_box(tmpdir):
    """
    Tests the absolute box size for a system with ligands far from
    the box
    """
    from vmd import atomsel, molecule
    from dabble import DabbleBuilder

    # Build the system
    p = str(tmpdir.mkdir("absolute_box"))
    b = DabbleBuilder(solute_filename=os.path.join(dir, "dor_ligs.mae"),
                      output_filename=os.path.join(p, "test.mae"),
                      user_x=75., user_y=75., user_z=115.,
                      overwrite=True, tmp_dir=p,
                      exclude_sel="same fragment as resname FRAG")
    b.write()

    # Load the built system
    m2 = molecule.load("mae", os.path.join(p, "test.mae"))
    molecule.set_top(m2)

    # Check all the ligands are there
    assert len(set(atomsel("resname FRAG").residue)) == 3

#==============================================================================

#def test_absolute_box_noexclude(tmpdir):
#    """
#    Tests the absolute box size for a system with ligands far from
#    the box, without using an exclude selection
#    """
#    from vmd import atomsel, molecule
#
#    from dabble import DabbleBuilder
#    p = str(tmpdir.mkdir("absolute_box"))
#
#    # Building the system should raise a valueerror in sanity check
#    # as resids are duplicated in protein chain
#    with pytest.raises(ValueError):
#        b = DabbleBuilder(solute_filename=os.path.join(dir, "dor_ligs.mae"),
#                          output_filename=os.path.join(p, "test.mae"),
#                          user_x=75., user_y=75., user_z=115.,
#                          overwrite=True, tmp_dir=p)
#        b.write()
#
#==============================================================================
