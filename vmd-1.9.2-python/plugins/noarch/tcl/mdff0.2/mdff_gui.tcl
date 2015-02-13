#package requre mdff
package require ssrestraints
package require namdrun

package provide mdff_gui 0.1
namespace eval MDFFGUI:: {
  namespace export mdffgui
  namespace eval gui { 
    variable w 
    variable CurrentPDBFile
    variable CurrentPSFFile
    variable DensityID
    variable CButtonChi
    variable MapID
    variable CButtonCis
    variable CCResInput
    variable CCSelInput
    variable CurrentDXFile
    variable Extrabonds 
    variable CurrentRefs
    
    variable CurrentCrystPDB
    
    variable CurrentDirLabel
    variable MolID
    variable env
    
    variable CCStep
    variable CCThreshInput
    variable UseCCThresh

    variable IMDSubmit
    variable IMDConnect
    variable IMDPause
    variable IMDFinish
    variable IMDStatusLabel
    variable IMDProcInput   
    variable IMDCCPlot
    variable Notebook
    variable justChanged
    variable GridpdbGen
    variable FixedpdbGen
    variable MDFFSetup
    variable GenerateXBonds
    
    variable CrystOption
    variable CurrentRefsEntry
    variable RefStepsInput
    variable CurrentCrystEntry
    variable GetRefs
    variable GetCryst
    variable CButtonMask
    variable CButtonBFS
    variable MaskResInput
    variable MaskCutoffInput
    
    variable IMDInpOrList
    variable IMDHost
    variable namdbin
    variable CurrentServer 
    variable ServerList
    variable IMDHostInput
    variable IMDServer 
    variable CButtonIMDWait
    variable CButtonIMDIgnore
    variable IMDPortInput
    variable IMDFreqInput
    variable IMDKeepInput
    variable HLF
    variable CCx
    variable CCy
    variable ShowBTNPadX
  }
  namespace eval settings:: {
    namespace export settings
    variable CurrentPDBPath ""
    variable CurrentPSFPath ""
    variable CurrentDXPath ""
    variable CurrentRefsPath ""
    variable CurrentDir [pwd]
    variable CurrentCrystPath ""
    variable Temperature 300
    variable FTemperature 300
    variable Minsteps 200
    variable Numsteps 50000
    variable GScale 0.3
    variable GridPDBSelText "protein and noh"
    variable FixedPDBSelText "none"
    variable FixedColumn "occupancy"
    variable SimulationName "mdff"
    variable SimulationStep "1"
    variable SSRestraints 0
    variable ChiralityRestraints 0
    variable CispeptideRestraints 0
    variable xMDFF 0
    variable BFS 0
    variable Mask 0
    variable MaskRes 5
    variable MaskCutoff 5 
    variable RefSteps 20000
    variable PBCorGBIS ""
    variable IMD 0
    variable IMDPort 2000
    variable IMDFreq 1
    variable IMDWait 0
    variable IMDIgnore 0
    variable IMDProcs ""
    variable IMDKeep 0
    variable CrossCorrelation 0
    variable CCSel "protein and noh"
    variable CCRes ""
    variable UseCCThresh 0
    variable CCThresh ""
    variable GridforceLite 0
    variable ParameterList [list [file join $env(CHARMMPARDIR) par_all27_prot_lipid_na.inp]]
  }

  namespace eval servers {} ;# create the "servers" namespace
}

proc MDFFGUI::gui::mdffgui {} {
  variable w 
  variable CurrentPDBFile
  variable CurrentPSFFile
  variable DensityID
  variable CButtonChi
  variable MapID
  variable CButtonCis
  variable CCResInput
  variable CCSelInput
  variable CurrentDXFile
  variable Extrabonds 
  variable CurrentRefs
  variable NumFrames 0
  
  variable IMDProcInput   
  variable CurrentCrystPDB
  
  variable CurrentDirLabel
  variable MolID
  variable env
  
  variable CCStep
  variable CCThreshInput
  variable UseCCThresh

  variable IMDSubmit
  variable IMDConnect
  variable IMDPause
  variable IMDFinish
  variable IMDStatusLabel
  variable IMDCCPlot
  variable Notebook
  variable justChanged
  variable GridpdbGen
  variable FixedpdbGen
  variable MDFFSetup
  variable GenerateXBonds
  
  variable CrystOption
  variable CurrentRefsEntry
  variable RefStepsInput
  variable CurrentCrystEntry
  variable GetRefs
  variable GetCryst
  variable CButtonMask
  variable CButtonBFS
  variable MaskResInput
  variable MaskCutoffInput
  
  variable IMDInpOrList
  variable IMDHost
  variable namdbin
  variable CurrentServer 
  variable ServerList
  variable IMDHostInput
  variable IMDServer 
  variable CButtonIMDWait
  variable CButtonIMDIgnore
  variable IMDPortInput
  variable IMDFreqInput
  variable IMDKeepInput
  variable HLF
  variable ShowBTNPadX
  
  set downPoint \u25BC
  set rightPoint \u25B6
  
  variable CCx
  variable CCy
  
  set Extrabonds {}


  set CrystOption ""

  set IMDHost "localhost"
  set namdbin "namd2"

  set MolID ""
  set MapID ""  

  set CCx ""
  set CCy ""
  set CCStep ""
  set justChanged 0
  
  set ShowBTNPadX 0
 
  # setup the theme depending on what is available
  set themeList [ttk::style theme names]
  if { [lsearch -exact $themeList "aqua"] != -1 } {
      ttk::style theme use aqua
      set ShowBTNPadX 18
  } elseif { [lsearch -exact $themeList "clam"] != -1 } {
      ttk::style theme use clam
  } elseif { [lsearch -exact $themeList "classic"] != -1 } {
      ttk::style theme use classic
  } else {
      ttk::style theme use default
  }
  MDFFGUI::gui::add_server "default localhost" {
    jobtype  local
    namdbin  {namd2 +p%d}
    maxprocs 1
    numprocs 1
    timeout  20

  }
 # set CurrentServer ""
  #source [file join $env(MDFFDIR) settings.tcl]
  set ServerList [list]
  foreach s [lsort [info vars ::MDFFGUI::servers::*]] {
    lappend ServerList [namespace tail $s]
  }
  set CurrentServer [lindex $ServerList 0]

  set w [toplevel ".mdffgui"]
  wm title $w "MDFF GUI"
  #wm resizable $w 1 1
  grid columnconfigure $w 0 -weight 1
  grid rowconfigure $w 0 -weight 1
  
  set HLF [ttk::frame $w.hlf]
  grid $w.hlf -column 0 -row 0 -sticky nsew
  # allow hlf to resize with window
  grid columnconfigure $w.hlf 0 -weight 1
  grid rowconfigure $w.hlf 1 -weight 1
 
  set FileMenu [ttk::menubutton $w.hlf.filemenu -text "File" -menu $w.hlf.filemenu.save -width 5]
  set FileMenuSave [menu $w.hlf.filemenu.save -tearoff no] 
  $FileMenuSave add command -label "Save Settings..." -command {MDFFGUI::gui::save_settings}
  $FileMenuSave add command -label "Load Settings..." -command {MDFFGUI::gui::load_settings}
   
  
  set Notebook [ttk::notebook $w.hlf.n]
  set NBTab1 [ttk::frame $w.hlf.n.f1]; # first page, which would get widgets gridded into it 
  $Notebook add $NBTab1 -text "MDFF Setup"
  ttk::notebook::enableTraversal $Notebook
  
  grid columnconfigure $w.hlf.n 0 -weight 1
  grid rowconfigure $w.hlf.n 0 -weight 1
  
  grid columnconfigure $w.hlf.n.f1 0 -weight 1
  grid rowconfigure $w.hlf.n.f1 0 -weight 1
  #MDFF Files
  set MFrame [ttk::frame $w.hlf.n.f1.main]
  grid columnconfigure $MFrame 0 -weight 1
  grid rowconfigure $MFrame 6 -weight 1

  set FileFrame [ttk::labelframe $w.hlf.n.f1.main.fileframe -labelanchor nw]
  grid columnconfigure $FileFrame 1 -weight 1
  #grid rowconfigure $FileFrame 0 -weight 1
  #set BFrame [frame $w.hlf.n.f1.main.bframe -relief raised -bd 2]
  set CurrentPDBFile [ttk::label $w.hlf.n.f1.main.fileframe.pdblabel -text "PDB File:"]
  set CurrentPSFFile [ttk::label $w.hlf.n.f1.main.fileframe.psflabel -text "PSF File:"]
  set CurrentPDBEntry [ttk::entry $w.hlf.n.f1.main.fileframe.pdbentry -textvariable MDFFGUI::settings::CurrentPDBPath -width 40]
  set CurrentPSFEntry [ttk::entry $w.hlf.n.f1.main.fileframe.psfentry -textvariable MDFFGUI::settings::CurrentPSFPath -width 40]
  set SelectPDB [ttk::button $w.hlf.n.f1.main.fileframe.button5 -text "Browse" -command {MDFFGUI::gui::get_pdb} ]
  set SelectPSF [ttk::button $w.hlf.n.f1.main.fileframe.button6 -text "Browse" -command {MDFFGUI::gui::get_psf} ]
  set LoadFileSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.seperator  -orient horizontal]   
  set LoadFiles [ttk::button $w.hlf.n.f1.main.fileframe.button7 -text "Load PSF/PDB" -command {MDFFGUI::gui::load_struct}]

  set CurrentDirLabel [ttk::label $w.hlf.n.f1.main.fileframe.dirlabel -text "Working Directory:"]
  set CurrentDirEntry [ttk::entry $w.hlf.n.f1.main.fileframe.direntry -textvariable MDFFGUI::settings::CurrentDir -width 40]
  set SelectDir [ttk::button $w.hlf.n.f1.main.fileframe.selectdir -text "Browse" -command {MDFFGUI::gui::get_dir} ] 
  set CurrentDirSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.cdseperator  -orient horizontal]   
  
  set RestraintFrame [ttk::frame $w.hlf.n.f1.main.fileframe.restraints]
  set RestraintLabel [ttk::label $w.hlf.n.f1.main.fileframe.restraints.label -text "Restraints: "]
  #grid columnconfigure $FileFrame 1 -weight 1
  set CButtonSecondary [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.ss -text "Secondary Structure Restraints" -variable MDFFGUI::settings::SSRestraints]
  set CButtonChirality [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.chi -text "Chirality Restraints" -variable MDFFGUI::settings::ChiralityRestraints]
  set CButtonCispeptide [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.cis -text "Cispeptide Restraints" -variable MDFFGUI::settings::CispeptideRestraints]
  set GenerateXBonds [ttk::button $w.hlf.n.f1.main.fileframe.xbonds -text "Generate Restraints" -command {MDFFGUI::gui::generate_xbonds} -state disabled]
  set XBondsSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.xbondsseperator  -orient horizontal]   
  
  set GridpdbLabel [ttk::label $w.hlf.n.f1.main.fileframe.label2 -text "Gridpdb selection text:"]
  set GridpdbInput [ttk::entry $w.hlf.n.f1.main.fileframe.entry2 -textvariable MDFFGUI::settings::GridPDBSelText -width 40]
  set GridpdbGen [ttk::button $w.hlf.n.f1.main.fileframe.gridpdb -text "Generate GridPDB" -command {MDFFGUI::gui::make_gridpdb} -state disabled]
  set GridpdbSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.gridpdbseperator  -orient horizontal]   

  set FixedpdbLabel [ttk::label $w.hlf.n.f1.main.fileframe.fixedlabel -text "Fixed PDB selection text:"]
  set FixedpdbInput [ttk::entry $w.hlf.n.f1.main.fileframe.fixedentry -textvariable MDFFGUI::settings::FixedPDBSelText -width 40]
  set FixedColumnLabel [ttk::label $w.hlf.n.f1.main.fileframe.fixedclabel -text "Fixed PDB Column:"]
  set FixedColumnInput [ttk::entry $w.hlf.n.f1.main.fileframe.fixedcentry -textvariable MDFFGUI::settings::FixedColumn -width 40]
  set FixedpdbGen [ttk::button $w.hlf.n.f1.main.fileframe.fixedpdb -text "Generate Fixed PDB" -command {MDFFGUI::gui::make_fixedpdb} -state disabled]
  set FixedpdbSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.fixedpdbseperator  -orient horizontal]   
 
  set ParamListLabel [ttk::label $w.hlf.n.f1.main.fileframe.paramlistlabel -text "Parameter Files: "]
  set ParamListScroll [ttk::scrollbar $w.hlf.n.f1.main.fileframe.paramlistscroll -command "$w.hlf.n.f1.main.fileframe.paramlistbox yview" -orient vertical]
  set ParamListBox [listbox $w.hlf.n.f1.main.fileframe.paramlistbox -height 5 -yscrollcommand "$ParamListScroll set" -listvariable MDFFGUI::settings::ParameterList]
  set ParamListAdd [ttk::button $w.hlf.n.f1.main.fileframe.paramlistadd -text "Add" -command {MDFFGUI::gui::add_paramfile}]
  set ParamListRemove [ttk::button $w.hlf.n.f1.main.fileframe.paramlistremove -text "Remove" -command {
   foreach ind [.mdffgui.hlf.n.f1.main.fileframe.paramlistbox curselection] {set MDFFGUI::settings::ParameterList [lreplace $MDFFGUI::settings::ParameterList $ind $ind]} 
  }]
  
  set CurrentDXFile [ttk::label $w.hlf.n.f1.main.fileframe.dxlabel -text "Density Map:"]
  set CurrentDXEntry [ttk::entry $w.hlf.n.f1.main.fileframe.dxentry -textvariable MDFFGUI::settings::CurrentDXPath -width 40]
  set SelectDX [ttk::button $w.hlf.n.f1.main.fileframe.dx -text "Browse" -command {MDFFGUI::gui::get_density} ]
  set MakeDX [ttk::button $w.hlf.n.f1.main.fileframe.griddx -text "Generate MDFF Potential" -command {MDFFGUI::gui::make_griddx} ]
  set DXSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.dxseperator  -orient horizontal]   
 
  
  set NameLabel [ttk::label $w.hlf.n.f1.main.fileframe.simnamelabel -text "Simulation Output Name:"]
  set NameInput [ttk::entry $w.hlf.n.f1.main.fileframe.entry3 -textvariable MDFFGUI::settings::SimulationName -width 40]
  
  set StepLabel [ttk::label $w.hlf.n.f1.main.fileframe.simsteplabel -text "Simulation Output Step:"]
  set StepInput [ttk::entry $w.hlf.n.f1.main.fileframe.stepentry -textvariable MDFFGUI::settings::SimulationStep -width 40]
  
  set SetupSeparator [ttk::separator $w.hlf.n.f1.main.setupseperator  -orient horizontal]
#there is code in place to allow this to start disabled and only open once a psf/pdb and map are
  #specified. however, this doesn't account for xMDFF where the user just puts in a name and
  #doesn't click browse, so this button is normal now by default
  set MDFFSetup [ttk::button $w.hlf.n.f1.main.setup -text "Generate NAMD files" -command {MDFFGUI::gui::mdff_setup} -state normal]
 
  set ShowMDFF [ttk::label $w.hlf.n.f1.main.showmdff -text "$rightPoint MDFF Files..." -anchor w]
  set HideMDFF [ttk::label $w.hlf.n.f1.main.fileframe.hidemdff -text "$downPoint MDFF Files" -anchor w]
  $FileFrame configure -labelwidget $HideMDFF
  bind $HideMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.fileframe
      grid .mdffgui.hlf.n.f1.main.showmdff
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showmdff 
      grid .mdffgui.hlf.n.f1.main.fileframe -row 0 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.fileframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  grid $HLF -row 0 -column 0 -sticky nswe
  grid $FileMenu -row 0 -column 0 -sticky nsw
  grid $Notebook -row 1 -column 0 -sticky nsew
  grid $MFrame -row 0 -column 0 -sticky nswe
  grid $ShowMDFF -row 0 -column 0 -sticky nswe -pady 5 -padx $ShowBTNPadX
 # grid $HideMDFF -row 0 -column 0 -sticky nsew
  #grid $BFrame -row 1 -column 0 -columnspan 5 -sticky nsew
 
  grid $CurrentDirLabel -row 1 -column 0 -columnspan 2 -sticky nsew
  grid $CurrentDirEntry -row 1 -column 1 -columnspan 2 -sticky nsew
  grid $SelectDir -row 1 -column 3 -sticky nsew
  #grid $CurrentDirSeparator -row 2 -column 0 -columnspan 5 -sticky ew
 
  grid $CurrentPSFFile -row 3 -column 0 -sticky nswe
  grid $CurrentPSFEntry -row 3 -column 1 -columnspan 2 -sticky nsew
  grid $CurrentPDBFile -row 4 -column 0 -sticky nswe
  grid $CurrentPDBEntry -row 4 -column 1 -columnspan 2 -sticky nsew
  grid $SelectPSF -row 3 -column 3 -sticky nswe
  grid $SelectPDB -row 4 -column 3 -sticky nswe
  #grid $LoadFileSeparator -row 6 -column 0 -columnspan 5 -sticky ew
  grid $LoadFiles -row 5 -column 0 -columnspan 5 -sticky nswe -padx 8 -pady 8

  
  grid $CurrentDXFile -row 7 -column 0 -sticky nsew
  grid $CurrentDXEntry -row 7 -column 1 -columnspan 2 -sticky nsew 
  grid $SelectDX -row 7 -column 3 -sticky nsew
  #grid $MakeDX -row 8 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
  #grid $DXSeparator -row 9 -column 0 -columnspan 5 -sticky ew

  grid $GridpdbLabel -row 10 -column 0 -sticky nsew
  grid $GridpdbInput -row 10 -column 1 -columnspan 2 -sticky nsew
 # grid $GridpdbGen -row 11 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
 # grid $GridpdbSeparator -row 12 -column 0 -columnspan 5 -sticky ew
  

  grid $FixedpdbLabel -row 13 -column 0 -sticky nsew
  grid $FixedpdbInput -row 13 -column 1 -columnspan 2 -sticky nsew
  grid $FixedColumnLabel -row 14 -column 0 -sticky nsew
  grid $FixedColumnInput -row 14 -column 1 -columnspan 2 -sticky nsew
  #grid $FixedpdbGen -row 15 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
 # grid $FixedpdbSeparator -row 16 -column 0 -columnspan 5 -sticky ew


  grid $NameLabel -row 17 -column 0 -sticky nsew
  grid $NameInput -row 17 -column 1 -columnspan 2 -sticky nsew
  
  grid $StepLabel -row 18 -column 0 -sticky nsew
  grid $StepInput -row 18 -column 1 -columnspan 2 -sticky nsew
  #grid $XBondsSeparator -row 18 -column 0 -columnspan 5 -sticky ew
 
  grid $RestraintFrame -row 19 -column 0 -sticky nsw -columnspan 4
  grid $RestraintLabel -row 0 -column 0 -sticky nsew
  grid $CButtonSecondary -row 0 -column 1 -sticky nsew
  grid $CButtonChirality -row 0 -column 2 -sticky nsew
  grid $CButtonCispeptide -row 0 -column 3 -sticky nsew
  #grid $GenerateXBonds -row 18 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
  
  grid $ParamListLabel -row 20 -column 0 -sticky nsew
  grid $ParamListBox -row 21 -column 0 -sticky nsew -columnspan 2 -rowspan 2
  grid $ParamListAdd -row 21 -column 3 -sticky nsew
  grid $ParamListRemove -row 22 -column 3 -sticky nsew
  grid $ParamListScroll -row 21 -column 2 -sticky nsew -rowspan 2
  
  grid $SetupSeparator -row 4 -column 0 -columnspan 5 -sticky ew
  grid $MDFFSetup -row 5 -column 0 -sticky nsew -padx 8 -pady 8
 
 
  #Parameters

 # set NBTab2 [ttk::frame $w.hlf.n.f2]; # second page
 # $Notebook add $NBTab2 -text "Parameters" 
  
  set ParamFrame [ttk::labelframe $w.hlf.n.f1.main.paramframe -labelanchor nw]
  grid columnconfigure $ParamFrame 1 -weight 1

  #grid columnconfigure $w.n.f2.main2 0 -weight 1
  set TempLabel [ttk::label $w.hlf.n.f1.main.paramframe.templabel -text "Temperature (K):"]
  set TempInput [ttk::entry $w.hlf.n.f1.main.paramframe.tempentry -textvariable MDFFGUI::settings::Temperature]
  
  set FTempLabel [ttk::label $w.hlf.n.f1.main.paramframe.ftemplabel -text "Final Temperature (K):"]
  set FTempInput [ttk::entry $w.hlf.n.f1.main.paramframe.ftempentry -textvariable MDFFGUI::settings::FTemperature]
  
  set MinStepLabel [ttk::label $w.hlf.n.f1.main.paramframe.minsteplabel -text "Minimization Steps:"]
  set MinStepInput [ttk::entry $w.hlf.n.f1.main.paramframe.minstepentry -textvariable MDFFGUI::settings::Minsteps]
  
  set NumStepLabel [ttk::label $w.hlf.n.f1.main.paramframe.numsteplabel -text "Time steps:"]
  set NumStepInput [ttk::entry $w.hlf.n.f1.main.paramframe.numstepentry -textvariable MDFFGUI::settings::Numsteps]
 
  set GScaleLabel [ttk::label $w.hlf.n.f1.main.paramframe.gscalelabel -text "Grid Scaling Factor:"]
  set GScaleInput [ttk::entry $w.hlf.n.f1.main.paramframe.gscaleentry -textvariable MDFFGUI::settings::GScale]

  set CButtonLite [ttk::checkbutton $w.hlf.n.f1.main.paramframe.gridforcelite -text "gridforcelite" -variable MDFFGUI::settings::GridforceLite]

  set EnvironFrame [ttk::frame $w.hlf.n.f1.main.paramframe.envframe]
  set EnvironLabel [ttk::label $w.hlf.n.f1.main.paramframe.envframe.label -text "System Environment: "]
  set PBC [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.pbc -variable MDFFGUI::settings::PBCorGBIS -value "-pbc" -text "Periodic Boundary Conditions"]
  set GBIS [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.gbis -variable MDFFGUI::settings::PBCorGBIS -value "-gbis" -text "Implicit Solvent"]
  
  set VAC [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.vac -variable MDFFGUI::settings::PBCorGBIS -value "" -text "Vacuum"]
  
  set ShowParams [ttk::label $w.hlf.n.f1.main.showparams -text "$rightPoint Simulation Parameters..." -anchor w]
  set HideParams [ttk::label $w.hlf.n.f1.main.paramframe.hideparams -text "$downPoint Simulation Parameters" -anchor w]
  $ParamFrame configure -labelwidget $HideParams
  bind $HideParams <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.paramframe
      grid .mdffgui.hlf.n.f1.main.showparams
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowParams <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showparams
      grid .mdffgui.hlf.n.f1.main.paramframe -row 1 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.paramframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
   
  #pack $Notebook
  #grid $MFrameP -row 0 -column 0

  grid $ShowParams -row 1 -column 0 -sticky w -pady 5 -padx $ShowBTNPadX
  grid $TempLabel -row 0 -column 0 -sticky nw
  grid $TempInput -row 0 -column 1 -sticky nsew
  
  grid $FTempLabel -row 1 -column 0 -sticky nw
  grid $FTempInput -row 1 -column 1 -sticky nsew

  grid $MinStepLabel -row 2 -column 0 -sticky nw
  grid $MinStepInput -row 2 -column 1 -sticky nsew
  
  grid $NumStepLabel -row 3 -column 0 -sticky nw
  grid $NumStepInput -row 3 -column 1 -sticky nsew

  grid $GScaleLabel -row 4 -column 0 -sticky nw
  grid $GScaleInput -row 4 -column 1 -sticky nsew
  
  grid $CButtonLite -row 5 -column 0 -sticky nw

  grid $EnvironFrame -row 6 -column 0 -sticky nsew -columnspan 4
  grid $EnvironLabel -row 0 -column 0 -sticky nsew
  grid $VAC -row 0 -column 1 -sticky nsew
  grid $PBC -row 0 -column 2 -sticky nsew
  grid $GBIS -row 0 -column 3 -sticky nsew

  #xMDFF 
  #set NBTab3 [ttk::frame $w.hlf.n.f3];
  #$Notebook add $NBTab3 -text "xMDFF" 


  set xMDFFFrame [ttk::labelframe $w.hlf.n.f1.main.xmdffframe ]
  grid columnconfigure $xMDFFFrame 1 -weight 1
  
  set CButtonxMDFF [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.xmdff -text "xMDFF" -variable MDFFGUI::settings::xMDFF -command {MDFFGUI::gui::xmdff_toggle}]
  
  set CurrentRefs [ttk::label $w.hlf.n.f1.main.xmdffframe.refs -text "Reflection Data File:"] 
  set CurrentRefsEntry [ttk::entry $w.hlf.n.f1.main.xmdffframe.refsentry -textvariable MDFFGUI::settings::CurrentRefsPath -width 40 -state disabled]
  set GetRefs [ttk::button $w.hlf.n.f1.main.xmdffframe.getrefs -text "Browse" -command {MDFFGUI::gui::get_refs} -state disabled ]
  
  set RefStepsLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.refstepslabel -text "Refinement Steps:"]
  set RefStepsInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.refstepsentry -textvariable MDFFGUI::settings::RefSteps -width 40 -state disabled]
  
  set CurrentCrystPDB [ttk::label $w.hlf.n.f1.main.xmdffframe.crystpdb -text "Symmetry PDB (optional):"] 
  set CurrentCrystEntry [ttk::entry $w.hlf.n.f1.main.xmdffframe.crystentry -textvariable MDFFGUI::settings::CurrentCrystPath -width 40 -state disabled]
  set GetCryst [ttk::button $w.hlf.n.f1.main.xmdffframe.getcrystpdb -text "Browse" -command {MDFFGUI::gui::get_cryst} -state disabled ]
  
  set CButtonMask [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.mask -text "Mask Density" -variable MDFFGUI::settings::Mask -state disabled -command {MDFFGUI::gui::mask_toggle}]
  
  set CButtonBFS [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.bfs -text "Calculate Beta Factors" -variable MDFFGUI::settings::BFS -state disabled]

  set MaskResLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.maskreslabel -text "Mask Resolution (A):"]
  set MaskResInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.maskresentry -textvariable MDFFGUI::settings::MaskRes -width 40 -state disabled]
  
  set MaskCutoffLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.maskcutofflabel -text "Mask Cutoff (A):"]
  set MaskCutoffInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.maskcutoffentry -textvariable MDFFGUI::settings::MaskCutoff -width 40 -state disabled]
  
  set ShowxMDFF [ttk::label $w.hlf.n.f1.main.showxmdff -text "$rightPoint xMDFF..." -anchor w]
  set HidexMDFF [ttk::label $w.hlf.n.f1.main.xmdffframe.hidexmdff -text "$downPoint xMDFF" -anchor w]
  $xMDFFFrame configure -labelwidget $HidexMDFF
  bind $HidexMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.xmdffframe
      grid .mdffgui.hlf.n.f1.main.showxmdff
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowxMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showxmdff
      grid .mdffgui.hlf.n.f1.main.xmdffframe -row 2 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.xmdffframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
  
  #grid $xMDFFFrame -row 0 -column 0
  grid $ShowxMDFF -row 2 -column 0 -sticky nsew -pady 5 -padx $ShowBTNPadX
  grid $CButtonxMDFF -row 0 -column 0 -sticky nw
  
  grid $CurrentRefs -row 1 -column 0 -sticky nw
  grid $CurrentRefsEntry -row 1 -column 1 -sticky nsew
  grid $GetRefs -row 1 -column 2 -sticky nsew

  grid $CurrentCrystPDB -row 2 -column 0 -sticky nw
  grid $CurrentCrystEntry -row 2 -column 1 -sticky nsew
  grid $GetCryst -row 2 -column 2 -sticky nsew
  
  grid $RefStepsLabel -row 3 -column 0 -sticky nw
  grid $RefStepsInput -row 3 -column 1 -sticky nsew

  grid $CButtonMask -row 4 -column 0 -sticky nsew
  grid $MaskResLabel -row 5 -column 0 -sticky nw
  grid $MaskResInput -row 5 -column 1 -sticky nsew
 
  grid $MaskCutoffLabel -row 6 -column 0 -sticky nw
  grid $MaskCutoffInput -row 6 -column 1 -sticky nsew
  
  grid $CButtonBFS -row 7 -column 0 -sticky nsew

  
  #IMDFF
  
  set IMDFrame [ttk::labelframe $w.hlf.n.f1.main.imdframe]
  grid columnconfigure $IMDFrame 1 -weight 1
  
  set CButtonIMD [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdff -text "IMD" -variable MDFFGUI::settings::IMD -command {MDFFGUI::gui::imd_toggle}]

  set IMDPortLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdportlabel -text "IMD Port:"]
  set IMDPortInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdportentry -textvariable MDFFGUI::settings::IMDPort -state disabled]
  
  set IMDFreqLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdfreqlabel -text "IMD Frequency:"]
  set IMDFreqInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdfreqentry -textvariable MDFFGUI::settings::IMDFreq -state disabled]
  
  set IMDKeepLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdkeeplabel -text "IMD Keep Frames:"]
  set IMDKeepInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdkeepentry -textvariable MDFFGUI::settings::IMDKeep -state disabled]
  
  set CButtonIMDWait [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdffwait -text "Wait for IMD Connection" -variable MDFFGUI::settings::IMDWait -state disabled]
  
  set CButtonIMDIgnore [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdffignore -text "Ignore IMD forces" -variable MDFFGUI::settings::IMDIgnore -state disabled]


  #set IMDInputSwitch [radiobutton $w.hlf.n.f4.main.imdradioinput -variable IMDInpOrList -value "input" -text "Input IMD Hostname" -command {MDFFGUI:select_imd_host}]
  #set IMDHostLabel [label $w.hlf.n.f4.main.imdhostlabel -text "IMD Hostname"]
  #set IMDHostInput [entry $w.hlf.n.f4.main.imdhostentry -textvariable IMDHost -state disabled]
  #set IMDHostLabel [label $w.hlf.n.f4.main.imdhostlabel -text "OR"]
 # set IMDListSwitch [radiobutton $w.hlf.n.f4.main.imdradiolist -variable IMDInpOrList -value "list" -text "IMD Host from List" -command {MDFFGUI:select_imd_host}]
  set IMDHostLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdhostlabel -text "IMD Server:"]
  set IMDServer [ttk::combobox $w.hlf.n.f1.main.imdframe.imdserver -textvariable MDFFGUI::gui::CurrentServer -values $MDFFGUI::gui::ServerList -state disabled]

  set IMDProcLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdproclabel -text "Processors:"]
  set IMDProcInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdprocentry -textvariable MDFFGUI::settings::IMDProcs -state disabled]
 
   
  bind $IMDServer <<ComboboxSelected>> {MDFFGUI::gui::imd_server}
  
  set ShowIMD [ttk::label $w.hlf.n.f1.main.showimd -text "$rightPoint IMD Parameters..." -anchor w]
  set HideIMD [ttk::label $w.hlf.n.f1.main.imdframe.hideimd -text "$downPoint IMD Parameters" -anchor w]
  $IMDFrame configure -labelwidget $HideIMD
  bind $HideIMD <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.imdframe
      grid .mdffgui.hlf.n.f1.main.showimd
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowIMD <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showimd
      grid .mdffgui.hlf.n.f1.main.imdframe -row 3 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.imdframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }

  #grid $MFrameI -row 0 -column 0

  grid $ShowIMD -row 3 -column 0 -sticky nsew -pady 5 -padx $ShowBTNPadX
    
  grid $CButtonIMD -row 0 -column 0 -sticky nw
  
  grid $CButtonIMDWait -row 1 -column 0 -sticky nw
  
  grid $CButtonIMDIgnore -row 2 -column 0 -sticky nw

  grid $IMDPortLabel -row 3 -column 0 -sticky nw
  grid $IMDPortInput -row 3 -column 1 -sticky nsew

  grid $IMDFreqLabel -row 4 -column 0 -sticky nw
  grid $IMDFreqInput -row 4 -column 1 -sticky nsew
  
  grid $IMDKeepLabel -row 5 -column 0 -sticky nw
  grid $IMDKeepInput -row 5 -column 1 -sticky nsew

  grid $IMDHostLabel -row 6 -column 0 -sticky nw
#  grid $IMDListSwitch -row 5 -column 3
  grid $IMDServer -row 6 -column 1 -sticky nsew
 
  grid $IMDProcLabel -row 7 -column 0 -sticky nw
  grid $IMDProcInput -row 7 -column 1 -sticky nsew
  
  #IMD Connect
  set NBTab4 [ttk::frame $w.hlf.n.f4];
  $Notebook add $NBTab4 -text "IMDFF Connect" 
  grid columnconfigure $NBTab4 0 -weight 1
  grid rowconfigure $NBTab4 0 -weight 1
  set IMDConnectFrame [ttk::frame $w.hlf.n.f4.main]
  grid columnconfigure $IMDConnectFrame 0 -weight 1
  set IMDButtonFrame [ttk::frame $w.hlf.n.f4.main.buttons]
  grid columnconfigure $IMDButtonFrame {0 1 2 3} -weight 1
  
  set IMDSubmit [ttk::button $w.hlf.n.f4.main.buttons.imdsubmit -text "Submit" -command {MDFFGUI::gui::submit_job} -state disabled ]
  
  set IMDConnect [ttk::button $w.hlf.n.f4.main.buttons.imdconnect -text "Connect" -command {MDFFGUI::gui::imd_connect} -state disabled ]
  
  set IMDPause [ttk::button $w.hlf.n.f4.main.buttons.imdpause -text "Pause" -command {MDFFGUI::gui::imd_pause} -state disabled ]
  
  set IMDFinish [ttk::button $w.hlf.n.f4.main.buttons.imdfinish -text "Finish" -command {MDFFGUI::gui::imd_kill} -state disabled ]
    
  set IMDStatusLabel [ttk::label $w.hlf.n.f4.main.imdstatuslabel -text "IMD Status:"]
  
  set CCFrame [ttk::labelframe $w.hlf.n.f4.main.ccframe -labelanchor nw]
  grid columnconfigure $CCFrame 1 -weight 1
  grid rowconfigure $CCFrame 3 -weight 1
  set CButtonCC [ttk::checkbutton $w.hlf.n.f4.main.ccframe.cc -text "Calculate Real-Time Cross Correlation" -variable MDFFGUI::settings::CrossCorrelation -command {MDFFGUI::gui::cc_toggle}]
  
  
  set CCSelLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccsellabel -text "Selection:"]
  set CCSelInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccselentry -textvariable MDFFGUI::settings::CCSel -width 40 -state disabled]
  
  set CCResLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccreslabel -text "Map Resolution:"]
  set CCResInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccresentry -textvariable MDFFGUI::settings::CCRes -width 40 -state disabled]
  
  set CButtonCCThresh [ttk::checkbutton $w.hlf.n.f4.main.ccframe.ccthresh -text "Use Threshold" -variable MDFFGUI::settings::UseCCThresh -command {MDFFGUI::gui::thresh_toggle}]
  set CCThreshLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccthreshlabel -text "Threshold:"]
  set CCThreshInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccthreshentry -textvariable MDFFGUI::settings::CCThresh -width 40 -state disabled]
  
  set IMDCCPlotFrame [ttk::frame $w.hlf.n.f4.main.ccframe.plot]
  grid columnconfigure $IMDCCPlotFrame 0 -weight 1
  set IMDCCPlot [multiplot embed $w.hlf.n.f4.main.ccframe.plot -xsize 600 -ysize 400 -title "Real-Time Cross Correlation" -xlabel "Timestep" -ylabel "Cross Correlation" -lines -linewidth 2 -marker point -radius 2 -autoscale ]
  set imdccplot [$IMDCCPlot getpath]
  set ShowCC [ttk::label $w.hlf.n.f4.main.showcc -text "$rightPoint Cross Correlation Analysis..." -anchor w]
  set HideCC [ttk::label $w.hlf.n.f4.main.ccframe.hidecc -text "$downPoint Cross Correlation Analysis" -anchor w]
  $CCFrame configure -labelwidget $HideCC
  bind $HideCC <Button-1> {
      grid remove .mdffgui.hlf.n.f4.main.ccframe
      grid .mdffgui.hlf.n.f4.main.showcc
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowCC <Button-1> {
      grid remove .mdffgui.hlf.n.f4.main.showcc
      grid .mdffgui.hlf.n.f4.main.ccframe -row 2 -column 0 -columnspan 4 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f4.main.ccframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
 # grid $IMDInputSwitch -row 5 -column 0 
 # grid $IMDHostInput -row 5 -column 1

  grid $IMDConnectFrame -row 0 -column 0 -sticky nsew
  grid $IMDButtonFrame -row 0 -column 0 -sticky nsew
  grid $IMDSubmit -row 0 -column 0 -sticky nsew
  grid $IMDConnect -row 0 -column 1 -sticky nsew
  grid $IMDPause -row 0 -column 2 -sticky nsew
  grid $IMDFinish -row 0 -column 3 -sticky nsew

  grid $IMDStatusLabel -row 1 -column 0 -columnspan 4 -sticky nsw
  grid $ShowCC -row 2 -column 0 -sticky nsew -pady 5 -columnspan 4 -padx $ShowBTNPadX
 
  grid $CButtonCC -row 0 -column 0 -sticky nsw -columnspan 4
   
  grid $CCSelLabel -row 1 -column 0 -sticky nsw
  grid $CCSelInput -row 1 -column 1 -sticky nsew
  
  grid $CCResLabel -row 2 -column 0 -sticky nsw
  grid $CCResInput -row 2 -column 1 -sticky nsw
  
  grid $CButtonCCThresh -row 3 -column 0 -sticky nsw
  grid $CCThreshLabel -row 4 -column 0 -sticky nsw
  grid $CCThreshInput -row 4 -column 1 -sticky nsw

  grid $IMDCCPlotFrame -row 5 -column 0 -sticky nsew -columnspan 4
  grid $imdccplot -row 0 -column 0 -sticky nsew
  
  #Density Map Tools
  #set NBTab5 [ttk::frame $w.hlf.n.f5];
  #$Notebook add $NBTab5 -text "Map Tools" 
  
  #Basic MDFF analysis
  #set NBTab6 [ttk::frame $w.hlf.n.f6];
  #$Notebook add $NBTab6 -text "Analysis" 
  
  
  bind $Notebook <<NotebookTabChanged>> {MDFFGUI::gui::resizeToActiveTab}
}

proc MDFFGUI::gui::mask_toggle {} {
  variable MaskResInput
  variable MaskCutoffInput

  if {$MDFFGUI::settings::Mask} {
    $MaskResInput configure -state normal  
    $MaskCutoffInput configure -state normal  
  } else {
    $MaskResInput configure -state disabled 
    $MaskCutoffInput configure -state disabled  
  }


}

proc MDFFGUI::gui::add_paramfile {} {
  set pathname [tk_getOpenFile -defaultextension "CHARMM Topology Files" -filetypes {{inp {.inp}} {top {.top}} {rtf {.rtf}} {all {*}}}]

  set MDFFGUI::settings::ParameterList [lappend MDFFGUI::settings::ParameterList $pathname]

}

proc MDFFGUI::gui::imd_server {} {  
  variable CurrentServer
  upvar MDFFGUI::servers::$CurrentServer currentserver
  set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
}

proc MDFFGUI::gui::save_settings {} {
  
  set pathname [tk_getSaveFile -defaultextension ".tcl" -filetypes {{tcl {.tcl}} {all {*}}}]

  if {$pathname != "" } {
    set f [open $pathname "w"] 
    foreach setting [lsort [info vars ::MDFFGUI::settings::*]] {
      puts $f "set $setting [list [set $setting]]" 
    }
    foreach server [lsort [info vars ::MDFFGUI::servers::*]] {
      upvar $server currentserver
      puts $f "MDFFGUI::gui::add_server \"[namespace tail $server]\" \{"
      foreach e [array names $server] {
        puts $f "  $e [list $currentserver($e)]"  
      }
      puts $f "\}"
    }
    close $f
  }
}

proc MDFFGUI::gui::load_settings {} {  
  set pathname [tk_getOpenFile -defaultextension ".tcl" -filetypes {{tcl {.tcl}} {all {*}}}]

  if {$pathname != "" } {
    source $pathname
    if {$::MDFFGUI::settings::CrossCorrelation} {MDFFGUI::gui::cc_toggle}
    if {$::MDFFGUI::settings::xMDFF} {MDFFGUI::gui::xmdff_toggle}
    if {$::MDFFGUI::settings::Mask} {MDFFGUI::gui::mask_toggle}
    if {$::MDFFGUI::settings::IMD} {MDFFGUI::gui::imd_toggle}
    if {$::MDFFGUI::settings::UseCCThresh} {MDFFGUI::gui::thresh_toggle}
    #loads structure and density...nice for me testing, but does a user want this?
    if {$::MDFFGUI::settings::CurrentDXPath != ""} {
      set MDFFGUI::gui::MapID [mol new $MDFFGUI::settings::CurrentDXPath]
      set MDFFGUI::gui::DensityID $MDFFGUI::gui::MapID  
    }
    if {$::MDFFGUI::settings::CurrentPSFPath != "" && $::MDFFGUI::settings::CurrentPDBPath !=""} {MDFFGUI::gui::load_struct}
  }
}

proc MDFFGUI::gui::resizeToActiveTab {} {
  variable Notebook
  variable HLF

  update idletasks

#  set dimW [winfo reqwidth [$Notebook select]]
#  set dimMinW [winfo reqwidth $Notebook]

  set dimH [winfo reqheight [$Notebook select]]
  set dimH [expr $dimH + 75]
  
#  if {$dimW > $dimMinW} {
#    set useW $dimW 
#  } else {
#    set useW $dimMinW
#  }

  #other W code for variable widths, but that doesn't look as good on MACOS so go with
  #static width for now
  set useW 800
  
  wm geometry $MDFFGUI::gui::w [format "%ix%i" $useW $dimH] 

}

proc MDFFGUI::gui::xmdff_toggle {} {
  variable CurrentRefsEntry
  variable RefStepsInput
  variable CurrentCrystEntry
  variable GetRefs
  variable GetCryst
  variable CButtonMask
  variable CButtonBFS

  if {$MDFFGUI::settings::xMDFF} {
    $CurrentRefsEntry configure -state normal
    $RefStepsInput configure -state normal
    $CurrentCrystEntry configure -state normal
    $GetRefs configure -state normal
    $GetCryst configure -state normal
    $CButtonMask configure -state normal
    $CButtonBFS configure -state normal
  } else {
    $CurrentRefsEntry configure -state disabled
    $RefStepsInput configure -state disabled
    $CurrentCrystEntry configure -state disabled
    $GetRefs configure -state disabled
    $GetCryst configure -state disabled
    $CButtonMask configure -state disabled
    $CButtonBFS configure -state disabled
  }

}

proc MDFFGUI::gui::thresh_toggle {} {
  variable CCThreshInput

  if {$MDFFGUI::settings::UseCCThresh} {
    $CCThreshInput configure -state normal
  } else {
    $CCThreshInput configure -state disabled
  }

}


proc MDFFGUI::gui::cc_toggle {} {
  variable CCResInput
  variable CCSelInput

  if {$MDFFGUI::settings::CrossCorrelation} {
    $CCSelInput configure -state normal
    $CCResInput configure -state normal
  } else {
    $CCSelInput configure -state disabled
    $CCResInput configure -state disabled
  }

}

proc MDFFGUI::gui::imd_kill {} {
  imd kill
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: NAMD job terminated."
}

proc MDFFGUI::gui::imd_pause {} {

  imd pause toggle
  #$IMDStatusLabel configure -text "IMD Status: Paused."
}

proc MDFFGUI::gui::imd_toggle {} {
  variable CButtonIMDWait
  variable CButtonIMDIgnore
  variable IMDPortInput
  variable IMDFreqInput
  variable IMDKeepInput
  variable IMDServer

  if {$MDFFGUI::settings::IMD} {
    $CButtonIMDWait configure -state normal
    $CButtonIMDIgnore configure -state normal
    $IMDPortInput configure -state normal
    $IMDFreqInput configure -state normal
    $IMDKeepInput configure -state normal
    $IMDServer configure -state normal
    $MDFFGUI::gui::IMDSubmit configure -state normal
    $MDFFGUI::gui::IMDConnect configure -state normal
    $MDFFGUI::gui::IMDPause configure -state normal
    $MDFFGUI::gui::IMDFinish configure -state normal
    $MDFFGUI::gui::IMDProcInput configure -state normal
  } else {
    $CButtonIMDWait configure -state disabled
    $CButtonIMDIgnore configure -state disabled
    $IMDPortInput configure -state disabled
    $IMDFreqInput configure -state disabled
    $IMDKeepInput configure -state disabled
    $IMDServer configure -state disabled
    $MDFFGUI::gui::IMDSubmit configure -state disabled
    $MDFFGUI::gui::IMDConnect configure -state disabled
    $MDFFGUI::gui::IMDPause configure -state disabled
    $MDFFGUI::gui::IMDFinish configure -state disabled
    $MDFFGUI::gui::IMDProcInput configure -state disabled
  }
}

proc MDFFGUI::gui::select_imd_host {} {
  variable IMDInpOrList
  variable IMDHostInput
  variable IMDServer 
  
  if {$IMDInpOrList == "input"} {
    $IMDServer configure -state disabled
    $IMDHostInput configure -state normal
  } elseif {$IMDInpOrList == "list"} {
    $IMDHostInput configure -state disabled
    $IMDServer configure -state normal
  }

}
proc MDFFGUI::gui::mdff_setup {} {
  
  
  variable Extrabonds
  
  #variable MDFFGUI::settings::CurrentDir 
  variable MolID
 
  
  
  variable CrystOption
  
 
  #variable MDFFGUI::settings::CurrentDir 
   
  #file rename -force $MDFFGUI::settings::CurrentDXPath $MDFFGUI::settings::CurrentDir 
#  if {$Extrabonds != 0} {
#    foreach extrabond $Extrabonds {
#      file rename -force $extrabond $MDFFGUI::settings::CurrentDir 
#    }
#  }
#  file rename -force "gridpdb.pdb" $MDFFGUI::settings::CurrentDir 
#  file rename -force "griddx.dx" $MDFFGUI::settings::CurrentDir 
  
  if {$MDFFGUI::settings::CurrentPSFPath == ""} {
    tk_messageBox -message "Please provide a PSF file."
  } elseif {$MDFFGUI::settings::CurrentPDBPath == ""} {
    tk_messageBox -message "Please provide a PDB file."
  } elseif {!$MDFFGUI::settings::xMDFF && $MDFFGUI::settings::CurrentDXPath == ""} {
    tk_messageBox -message "Please provide a density map file."
  } elseif {$MolID == ""} {
    tk_messageBox -message "Please load the PSF and PDB using the button."
  } elseif {$MDFFGUI::settings::xMDFF && $MDFFGUI::settings::CurrentRefsPath == ""} {
    tk_messageBox -message "Please provide a reflections data file."
  } else {
    set refs ""
    set refsteps ""
    set xmdff ""
    set bfs ""
    set mask ""
    set maskres ""
    set maskcutoff ""

    
    file mkdir $MDFFGUI::settings::CurrentDir 
    
    if {[llength $Extrabonds] == 0} {set Extrabonds 0}
    if {$MDFFGUI::settings::xMDFF} {
      set xmdff "--xmdff"
      set refsteps "-refsteps $MDFFGUI::settings::RefSteps"
      #set refs "-refs [lindex [split $MDFFGUI::settings::CurrentRefsPath "/"] end]"
      set refs "-refs [file tail $MDFFGUI::settings::CurrentRefsPath]"
      set maskres "-mask_res $MDFFGUI::settings::MaskRes"
      set maskcutoff "-mask_cutoff $MDFFGUI::settings::MaskCutoff"
      file copy -force $MDFFGUI::settings::CurrentRefsPath $MDFFGUI::settings::CurrentDir 
      if {$CrystOption != ""} {
        file copy -force $MDFFGUI::settings::CurrentCrystPath $MDFFGUI::settings::CurrentDir 
      }
      if {$MDFFGUI::settings::BFS} {
        set bfs "--bfs"
      }
      if {$MDFFGUI::settings::Mask} {
       set mask "--mask"
      }
    }

    set fixedpdb ""
    set fixedcol ""

    if {$MDFFGUI::settings::FixedPDBSelText != "none"} {
      set fixedpdb "-fixpdb fixed.pdb"
      set fixedcol "-fixcol O"
#    file rename -force "fixed.pdb" MDFFGUI::$MDFFGUI::settings::CurrentDir 
    }
     
    set parlist ""
    foreach par $MDFFGUI::settings::ParameterList {
      file copy -force $par $MDFFGUI::settings::CurrentDir
      set parlist [lappend parlist [file tail $par]]
    }
    set imd ""
    set imdport ""
    set imdfreq ""
    set imdwait ""
    set imdignore ""
    if {$MDFFGUI::settings::IMD} {
      set imd "--imd"
      set imdport "-imdport $MDFFGUI::settings::IMDPort"
      set imdfreq "-imdfreq $MDFFGUI::settings::IMDFreq"
      if {$MDFFGUI::settings::IMDWait} {
        set imdwait "--imdwait"
      }
      if {$MDFFGUI::settings::IMDIgnore} {
        set imdignore "--imdignore"
      }
    }
    
    set lite ""
    if {$MDFFGUI::settings::GridforceLite} {
      set lite "--lite"
    }
    
    file copy -force $MDFFGUI::settings::CurrentPSFPath $MDFFGUI::settings::CurrentDir 
    file copy -force $MDFFGUI::settings::CurrentPDBPath $MDFFGUI::settings::CurrentDir 

    MDFFGUI::gui::make_gridpdb
    MDFFGUI::gui::make_fixedpdb
    if {!$MDFFGUI::settings::xMDFF} {
      MDFFGUI::gui::make_griddx
    }
    MDFFGUI::gui::generate_xbonds

    mdff setup -o $MDFFGUI::settings::SimulationName -dir $MDFFGUI::settings::CurrentDir -psf [file tail $MDFFGUI::settings::CurrentPSFPath] \
    -pdb [file tail $MDFFGUI::settings::CurrentPDBPath] -griddx "griddx.dx" \
    -gridpdb "gridpdb.pdb" -extrab $Extrabonds -gscale $MDFFGUI::settings::GScale -minsteps $MDFFGUI::settings::Minsteps -numsteps $MDFFGUI::settings::Numsteps -step 1 \
    -minsteps $MDFFGUI::settings::Minsteps -numsteps $MDFFGUI::settings::Numsteps -temp $MDFFGUI::settings::Temperature -ftemp $MDFFGUI::settings::FTemperature $MDFFGUI::settings::PBCorGBIS \
    $xmdff {*}$refs {*}$refsteps {*}$CrystOption $mask $bfs $maskres $maskcutoff $imd $imdport $imdfreq $imdwait $imdignore \
    {*}$fixedpdb {*}$fixedcol $lite -parfiles $parlist -step $MDFFGUI::settings::SimulationStep 
  }
}

proc MDFFGUI::gui::add_server {args} {
  variable ServerList
  variable CurrentServer
  variable IMDServer

  set servername "[lindex $args 0]"
  set opts [lindex $args 1]
  
  array set $servername $opts

#  lappend ServerList $servername
#   lappend ServerList [array get $servername]
 
  if {[info exists IMDServer] && [lsearch $ServerList $servername] == -1} {
    array set "::MDFFGUI::servers::$servername" $opts
    set ServerList [lappend ServerList $servername]
    $IMDServer configure -values $ServerList
    set CurrentServer [lindex $ServerList 0]
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
  } elseif {![info exists IMDServer]} {
    array set "::MDFFGUI::servers::$servername" $opts
    set ServerList [lappend ServerList $servername]
    set CurrentServer [lindex $ServerList 0]
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
  }
 #puts "adding server: $servername into $ServerList" 
 # set ServerList [list]
 # foreach s [lsort [info vars ::MDFFGUI::servers::*]] {
 #   lappend ServerList [namespace tail $s]
 # }
 # set CurrentServer [lindex $ServerList 0]
}

proc MDFFGUI::gui::get_dir {} {
  #variable MDFFGUI::settings::CurrentDir 
  variable CurrentDirLabel

  set pathname [tk_chooseDirectory]
  if {$pathname != ""} {
 #   set fname [lindex [split $pathname "/"] end]
 #   $CurrentDirLabel configure -text "Current Working Directory: $fname/"
    set MDFFGUI::settings::CurrentDir $pathname
  }
  
}

proc MDFFGUI::gui::submit_job {} {
  variable namdbin
  variable CurrentServer
  variable IMDHost
  variable CCStep
  variable CCx
  variable CCy
  #variable MDFFGUI::settings::CurrentDir 
  variable IMDInpOrList
  variable MolID
    
   mol top $MolID
   $MDFFGUI::gui::IMDCCPlot clear
   set CCx ""
   set CCy ""
  # $IMDCCPlot add $CCx $CCy
   $MDFFGUI::gui::IMDCCPlot replot
#    if {$IMDInpOrList == "list"}{   
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set namdbin "$currentserver(namdbin)"
    set exec_command "$namdbin $MDFFGUI::settings::SimulationName-step$MDFFGUI::settings::SimulationStep.namd"
#    } elseif {$IMDInpOrList == "input"}{
#      array set currentserver {
#        jobtype  local
#        namdbin  {namd2 +p%d}
#        maxprocs 12
#        numprocs 12
#        timeout  20
#      }
#    }
 # set currdir [exec pwd]
  set imdjobspec [NAMDRun::submitJob currentserver $exec_command $MDFFGUI::settings::CurrentDir $MDFFGUI::settings::IMDProcs automdff.log] 

  set IMDHost [lindex $imdjobspec 1]
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: NAMD job submitted."
}

proc MDFFGUI::gui::imd_update {args} {
  variable MolID
  variable IMDFreq
  variable justChanged
  variable CCx
  variable CCy
  variable CC
  variable CCStep
  variable NumFrames

  variable DensityID
  variable MapID
  #puts $xMDFF
  set ts [molinfo $MolID get timesteps]
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: Step $ts"
  set MapReps ""
  set MapColors ""
  set MapMats ""
  set CMMap ""
  set RMMap ""
  set SMMap ""
  set GMMap ""
  set CM ""
  set RM ""
  set SM ""
  set GM ""
  
  #xMDFF map loading
  if {$MDFFGUI::settings::xMDFF} {
    if {[expr ($ts-$MDFFGUI::settings::Minsteps-100) % $MDFFGUI::settings::RefSteps] == 0 && $justChanged != $ts || $ts == 1 && $justChanged != $ts} {
      if {$MapID != ""} {
        #set CMMap [molinfo $MapID get center_matrix]
        #set RMMap [molinfo $MapID get rotate_matrix]
        #set SMMap [molinfo $MapID get scale_matrix]
        #set GMMap [molinfo $MapID get global_matrix]
        for {set i 0} {$i < [molinfo $MapID get numreps] } {incr i} {
          set MapReps [lappend MapReps [lindex [molinfo $MapID get "{rep $i}"] 0]]
          set MapColors [lappend MapColors [lindex [molinfo $MapID get "{color $i}"] 0]]
          set MapMats [lappend MapMats [lindex [molinfo $MapID get "{material $i}"] 0]]
        }
        
        mol delete $MapID 
      }
      set CM [molinfo $MolID get center_matrix]
      set RM [molinfo $MolID get rotate_matrix]
      set SM [molinfo $MolID get scale_matrix]
      set GM [molinfo $MolID get global_matrix]
      set MapID [mol new [file join $MDFFGUI::settings::CurrentDir "xmdff_density.dx"]]
      mol top $MolID
      #molinfo $MolID set center_matrix $CM
      #molinfo $MolID set rotate_matrix $RM
      #molinfo $MolID set scale_matrix $SM
      #molinfo $MolID set global_matrix $GM
      #molinfo $MapID set center_matrix $CMMap
      #molinfo $MapID set rotate_matrix $RMMap
      #molinfo $MapID set scale_matrix $SMMap
      #molinfo $MapID set global_matrix $GMMap
      foreach mol [molinfo list] {
        molinfo $mol set center_matrix $CM 
        molinfo $mol set rotate_matrix $RM
        molinfo $mol set scale_matrix $SM
        molinfo $mol set global_matrix $GM
      }


      for {set i 0} {$i < [llength $MapReps]} {incr i} {
        if {$i > 0} { mol addrep $MapID }
        mol modstyle $i $MapID [lindex $MapReps $i]
        mol modcolor $i $MapID [lindex $MapColors $i]
        mol modmaterial $i $MapID [lindex $MapMats $i]
      }

      set justChanged $ts
    }
  }

  #CC plot
  if {$MDFFGUI::settings::CrossCorrelation} {
    MDFFGUI::gui::resizeToActiveTab
    if {$MDFFGUI::settings::CCRes != "" && $CCStep != $ts} {
      set sel [atomselect $MolID $MDFFGUI::settings::CCSel]
      if {$MDFFGUI::settings::xMDFF} { 
        set DensityMol $MapID 
      } else {
        set DensityMol $DensityID
      }
      if {!$MDFFGUI::settings::UseCCThresh} {
        set CCy [lappend CCy [mdffi cc $sel -res $MDFFGUI::settings::CCRes -mol $DensityMol]]
        set CCx [lappend CCx $ts]
      } elseif {[string is double $MDFFGUI::settings::CCThresh] || [string is integer $MDFFGUI::settings::CCThresh] } {
        set CCy [lappend CCy [mdffi cc $sel -res $MDFFGUI::settings::CCRes -mol $DensityMol -thresholddensity $MDFFGUI::settings::CCThresh]]
        set CCx [lappend CCx $ts]
      }
        $MDFFGUI::gui::IMDCCPlot clear
        $MDFFGUI::gui::IMDCCPlot add $CCx $CCy
        $MDFFGUI::gui::IMDCCPlot replot
      set CCStep $ts
    } else {
    }
  }

}

proc MDFFGUI::gui::imd_connect {} {
  variable MolID
  variable IMDCCPlot
  #variable vmd_timestep

  set attempt_delay 1500   ;# delay between attepmts to connect in ms
  set attempt_timeout 10000 ;# timeout in ms
   
   
  mol top $MolID
  #$IMDCCPlot clear
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: Trying to connect (waiting for NAMD)..."
  
  update idletasks
  for { set timecounter 0 } { $timecounter <= $attempt_timeout } {incr timecounter $attempt_delay} {
    if ![catch { imd connect $MDFFGUI::gui::IMDHost $MDFFGUI::settings::IMDPort }] {
      imd keep 0
      #setstate "connected"
      #trace add variable ::vmd_timestep($imdmol) write AutoMDFF::tracetimestep
      trace add variable ::vmd_timestep($MolID) write ::MDFFGUI::gui::imd_update   
      imd keep $MDFFGUI::settings::IMDKeep
      return
    }
    # else give NAMD more time
    after $attempt_delay
  }

}


proc MDFFGUI::gui::get_cryst {} {
  variable CurrentCrystPDB
  
  variable CrystOption

  set pathname [tk_getOpenFile -defaultextension ".pdb" -filetypes {{pdb {.pdb}} {all {*}}}]
  if {$pathname != ""} {
#    $CurrentCrystPDB configure -text "Symmetry PDB: $fname"
    set MDFFGUI::settings::CurrentCrystPath $pathname
    set CrystOption "-crystpdb [file tail $pathname]"
  }
}

proc MDFFGUI::gui::get_refs {} {
  variable CurrentRefs
  

  set pathname [tk_getOpenFile -defaultextension ".mtz, .cif" -filetypes {{cif {.cif}} {mtz {.mtz}} {all {*}}}]
  if {$pathname != ""} {
    #set fname [lindex [split $pathname "/"] end]
   # $CurrentRefs configure -text "Reflection Data File: $fname"
    set MDFFGUI::settings::CurrentRefsPath $pathname
  }
}

proc MDFFGUI::gui::make_gridpdb {} {
  file mkdir $MDFFGUI::settings::CurrentDir 
  mdff gridpdb -psf $MDFFGUI::settings::CurrentPSFPath -pdb $MDFFGUI::settings::CurrentPDBPath -o [file join $MDFFGUI::settings::CurrentDir "gridpdb.pdb"] -seltext $MDFFGUI::settings::GridPDBSelText
  mol delete top
}

proc MDFFGUI::gui::make_fixedpdb {} {
  #variable MDFFGUI::settings::CurrentDir 
  variable MolID
  
  file mkdir $MDFFGUI::settings::CurrentDir
  set selall [atomselect $MolID all] 
  set sel [atomselect $MolID $MDFFGUI::settings::FixedPDBSelText]
  $selall set $MDFFGUI::settings::FixedColumn 0
  $sel set $MDFFGUI::settings::FixedColumn 1
  $selall writepdb [file join $MDFFGUI::settings::CurrentDir "fixed.pdb"]
}

proc MDFFGUI::gui::make_griddx {} {
  
  #variable MDFFGUI::settings::CurrentDir 
  variable MDFFSetup
  variable MolID
  variable DensityID

  if {$MDFFGUI::settings::CurrentDXPath == ""} {
    tk_messageBox -message "Please provide a density map file."
  } else {
    file mkdir $MDFFGUI::settings::CurrentDir 
    mdff griddx -i $MDFFGUI::settings::CurrentDXPath -o [file join $MDFFGUI::settings::CurrentDir "griddx.dx"]
    $MDFFSetup configure -state normal
    set DensityID [mol new $MDFFGUI::settings::CurrentDXPath]
    mol inactive $DensityID
    mol off $DensityID
#    mol top $MolID
  }
}

proc MDFFGUI::gui::get_density {} {
  variable CurrentDXFile
  
  
  variable MDFFSetup 
  variable MolID
  variable MapID

  set pathname [tk_getOpenFile -defaultextension ".dx, .situs" -filetypes {{dx {.dx}} {situs {.situs}} {ccp4 {.ccp4}} {all {*}}}]
  if {$pathname != ""} {
#    $CurrentDXFile configure -text "Density Map: $fname"
    set MDFFGUI::settings::CurrentDXPath $pathname
    set MapID [mol new $MDFFGUI::settings::CurrentDXPath]
#    mol top $MolID
    if {$MDFFGUI::settings::CurrentPSFPath != "" && $MDFFGUI::settings::CurrentPDBPath != ""} {
      $MDFFSetup configure -state normal
    }
  }
}

proc MDFFGUI::gui::load_struct {} {
  variable MolID
  
  variable GridpdbGen
  variable FixedpdbGen
  
  variable MDFFSetup
  variable GenerateXBonds
  
  if {$MDFFGUI::settings::CurrentPSFPath == ""} {
    tk_messageBox -message "Please provide a PSF file."
  } elseif {$MDFFGUI::settings::CurrentPDBPath == ""} {
    tk_messageBox -message "Please provide a PDB file."
  } else {  
    set MolID [mol new $MDFFGUI::settings::CurrentPSFPath]
    mol addfile $MDFFGUI::settings::CurrentPDBPath
    
    if {$MDFFGUI::settings::CurrentDXPath != ""} {
      $MDFFSetup configure -state normal
    }
  }
#  $GridpdbGen configure -state normal
#  $FixedpdbGen configure -state normal
#  $GenerateXBonds configure -state normal
}

proc MDFFGUI::gui::generate_xbonds {} {
  
  variable MolID
  variable Extrabonds
  #variable MDFFGUI::settings::CurrentDir 
#  variable CButtonSS
#  variable CButtonCis
#  variable CButtonChi


#  variable CButtonSS
 # variable CButtonChi
 # variable CButtonCis
  
  set Extrabonds ""

  file mkdir $MDFFGUI::settings::CurrentDir 
#CHANGE TO JUST USE TOP MOL AND MAKE USE OF MOLINFO GET FILENAME
  if {$MDFFGUI::settings::SSRestraints} { 
    ssrestraints -psf $MDFFGUI::settings::CurrentPSFPath -pdb $MDFFGUI::settings::CurrentPDBPath -o [file join $MDFFGUI::settings::CurrentDir "ssrestraints.txt"] -hbonds
    lappend Extrabonds "ssrestraints.txt" 
    mol delete top 
  }
  if {$MDFFGUI::settings::CispeptideRestraints} { 
    mol top $MolID
    cispeptide restrain -o [file join $MDFFGUI::settings::CurrentDir "cispeptide.txt"]
    lappend Extrabonds "cispeptide.txt"
  }
  if {$MDFFGUI::settings::ChiralityRestraints} {
    mol top $MolID
    chirality restrain -o [file join $MDFFGUI::settings::CurrentDir "chirality.txt"]
    lappend Extrabonds "chirality.txt"
  }
  if {[llength $Extrabonds] == 0} {
    set Extrabonds 0
  }

}
#CHANGE TO JUST USE TOP MOL AND MAKE USE OF MOLINFO GET FILENAME
proc MDFFGUI::gui::get_pdb {} {
  variable CurrentPDBFile
 
  set pathname [tk_getOpenFile -defaultextension ".pdb" -filetypes {{pdb {.pdb}} {all {*}}}]
  if {$pathname != ""} {
#    $CurrentPDBFile configure -text "PDB File: $fname"
    set MDFFGUI::settings::CurrentPDBPath $pathname
  }
}

proc MDFFGUI::gui::get_psf {} {
  variable CurrentPSFFile
  

  set pathname [tk_getOpenFile -defaultextension ".psf" -filetypes {{psf {.psf}} {all {*}}}]
  if {$pathname != ""} {
 #   $CurrentPSFFile configure -text "PSF File: $fname"
    set MDFFGUI::settings::CurrentPSFPath $pathname
  }
}

proc mdffgui_tk {} {
  MDFFGUI::gui::mdffgui
  return $MDFFGUI::gui::w
}
