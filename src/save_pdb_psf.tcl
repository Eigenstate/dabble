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
set dabbledir $::env(DABBLEDIR)
topology ${dabbledir}/charmm_params/top_water_ions.rtf
topology ${dabbledir}/charmm_params/top_all36_cgenff.rtf
topology ${dabbledir}/charmm_params/top_all36_prot.rtf
topology ${dabbledir}/charmm_params/top_all36_lipid.rtf
puts "INFO: Done reading topology"

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
  # Residue substitutions
  # DONE IN PYTHON CODE NOW THAT PRODUCES INPUT PDB

  # Segment definition needs to come AFTER residue aliases!
  segment P${i} {
    pdb $protnam
  }
  
  #EXTRA_PROTEIN_COMMANDS
  coordpdb $protnam P${i}
  incr i
}

# Write the output files
writepsf x-plor cmap ${output}.psf
writepdb ${output}.pdb


