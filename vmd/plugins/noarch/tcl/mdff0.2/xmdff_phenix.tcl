package require mdff
set argv $env(VMDARG)
#psf file
set psf [lindex $argv 0]   
#OUTPUTNAME,INPUTNAME or PDB
set coord [lindex $argv 1]
#GRIDFILE, WARNING WILL OVERWRITE FILE
set maps [lindex $argv 2]
#diffraction data
set refs [lindex $argv 3]
#Do we want to use individual adp refinement?
set bfs [lindex $argv 4]
#Do we want to mask the map?
set mask [lindex $argv 5]
#pdb which never changes that contains a valid CRYST1 line for symmetry information
set crystpdb [lindex $argv 6]
#resolution of mask density
set mask_res [lindex $argv 7]
#cutoff for mask density
set mask_cutoff [lindex $argv 8]

file delete -force "mapinput.pdb"
file delete -force "[list $maps]"
file delete -force "mask.dx"
file delete -force "mapinput_2mFo-DFc_map.ccp4"
file delete -force "betainitial.pdb"
file delete -force "betainitial_refine_001.pdb"

mol new $psf waitfor all

if ([file exists "$coord.restart.coor"]) {
  mol addfile "$coord.restart.coor" waitfor all
} else {
  mol addfile "$coord" waitfor all
}

set sel [atomselect top "protein and noh"]

$sel frame last
$sel set occupancy 1
$sel writepdb "betainitial.pdb"

set frpdb [open "betainitial.pdb" "r"]
set spdb [read $frpdb]
close $frpdb
set fwpdb [open "betainitial.pdb" "w"]
regsub -all "HSD" $spdb "HIS" spdb
regsub -all "HSE" $spdb "HIS" spdb
regsub -all "URA" $spdb "  U" spdb
regsub -all "ADE" $spdb "  A" spdb
regsub -all "CYT" $spdb "  C" spdb
regsub -all "GUA" $spdb "  G" spdb
regsub -all "CYN" $spdb "CYS" spdb
regsub -all -line {^.*CRYST.*$} $spdb " " spdb
puts $fwpdb $spdb
close $fwpdb

if ($bfs) {
  puts "calculating beta factors..."
  #may not like symmetry format in cif/mtz.  may have to add a --symmetry="pdbwithCRYST1line.pdb" that never gets changed
  if {$crystpdb != 0} {
    exec phenix.refine "betainitial.pdb" "$refs" refinement.refine.strategy=individual_adp --symmetry=$crystpdb --overwrite
  } else {
    exec phenix.refine "betainitial.pdb" "$refs" refinement.refine.strategy=individual_adp --overwrite
  }
  exec cp "betainitial_refine_001.pdb" "mapinput.pdb"
} else {
  exec cp "betainitial.pdb" "mapinput.pdb"
}  

set frpdb [open "mapinput.pdb" "r"]
set spdb [read $frpdb]
close $frpdb
set fwpdb [open "mapinput.pdb" "w"]
regsub -all -line {^.*CRYST.*$} $spdb " " spdb
puts $fwpdb $spdb
close $fwpdb

for {set i 0} {$i < [llength $refs]} {incr i} {
  set map [lindex $maps [expr [llength $maps]-[llength $refs]+$i]]
  puts "computing density map..."
  exec phenix.maps "maps$i.params"

#in case of crashes due to R-frees
#  exec phenix.remove_free_from_map mapinput_map_coeffs.mtz 3p5Af_f.mtz
#  exec phenix.mtz2map mtz_file=map_coeffs_without_freer_set.mtz pdb_file=mapinput.pdb
#  mdff griddx -i map_coeffs_without_freer_set_2mFo-DFc.ccp4 -o "$map"

  file delete -force "xmdff_density$i.dx"
  mdff griddx -i mapinput_2mFo-DFc_map.ccp4 -o "$map"
  mdff griddx -i "$map" -o "xmdff_density$i.dx"

#masking begin
  if ([lindex $mask $i]) {
    puts "masking map..."
    volmap mask $sel -res [lindex $mask_res $i] -cutoff [lindex $mask_cutoff $i] -o "mask.dx" 
    volutil -mult "xmdff_density$i.dx" "mask.dx" -o "xmdff_density$i.dx"
    mdff griddx -i "xmdff_density$i.dx" -o "$map"
    #mdff griddx -i "$map" -o "xmdff_density.dx"
  }
#masking end
}

exit
