# Fix the package search path
#lappend auto_path [file join $env(VMDDIR) plugins noarch tcl]
#lappend auto_path [file join $env(VMDDIR) plugins [vmdinfo arch] tcl autopsf1.4]

# This is autopsf pkgIndex.tcl
set dir [file join $env(VMDDIR) plugins [vmdinfo arch] tcl psfgen1.6]
package ifneeded psfgen 1.6.2 [list load [file join $dir libpsfgen.so]]

#source [file join $env(VMDDIR) plugins [vmdinfo arch] tcl autopsf1.4 pkgIndex.tcl]
package require psfgen

#$inputpdb = [lindex $argv 0]
#$output = [lindex $argv 1]
set output "psfgen_output"
puts $output

# Load the parameters
resetpsf
topology /opt/vmd_src/plugins/readcharmmtop/top_all27_prot_lipid_na.inp

# First pull out the water and set the segment
set waterfiles [glob -directory /tmp psf_wat_*.pdb]
set i 0
foreach watnam $waterfiles {
  set mid mol new $watnam
  set water [atomselect top "water"]
  segment W${i} {
    auto none
    pdb $watnam
  }
  pdbalias atom HOH O OH2
  coordpdb $watnam W${i}
  mol delete $mid
  incr i
}

# TODO: Ions
set ionfile [glob -directory /tmp psf_ions_*.pdb]
set i 0
foreach ionnam $ionfile {
  set mid [mol new $ionfile]
  segment I${i} {
    pdb $ionnam
  }
  coordpdb $ionnam I${i}
  mol delete $mid
  incr i
}


# Now the lipids
set lipidfile [glob -directory /tmp psf_lipid_*.pdb]
set i 0
foreach lipnam $lipidfile {
  set mid [mol new $lipidfile]
  segment L${i} {
    pdb $lipnam
  }
  coordpdb $lipnam L${i}
  mol delete $mid
  incr i
}

# Now protein (assume one chain)
set protfile [glob -directory /tmp psf_protein_*.pdb]
set i 0
foreach protnam $protfile {
  segment P${i} {
    pdb $protnam
  }
  pdbalias residue HIS HSE
  pdbalias atom ILE CD1 CD
  coordpdb $protnam P${i}
  incr i
}

# Write the output files
writepsf ${output}.psf
writepdb ${output}.pdb


