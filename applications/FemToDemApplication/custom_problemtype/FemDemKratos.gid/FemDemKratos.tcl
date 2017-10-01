## GiD events --------------------------------------------------------------------------------------------------------------------------------------------------

proc InitGIDProject { dir } {
    
    # Initialize ProblemType Menu
    if { [GidUtils::IsTkDisabled] eq 0} {  
        GiDMenu::Create "FemDemKratos Application" PRE
        GiDMenu::InsertOption "FemDemKratos Application" [list "Parts"] 0 PRE "GidOpenConditions \"Parts\"" "" ""
        GiDMenu::InsertOption "FemDemKratos Application" [list "Dirichlet Constraints"] 1 PRE "GidOpenConditions \"Dirichlet_Constraints\"" "" ""
        GiDMenu::InsertOption "FemDemKratos Application" [list "Loads"] 2 PRE "GidOpenConditions \"Loads\"" "" ""
        GiDMenu::InsertOption "FemDemKratos Application" [list "Project Parameters"] 3 PRE "GidOpenProblemData" "" ""
        GiDMenu::UpdateMenus
    }

    # Save ProblemTypePath
    #set ::FemDemKratos::ProblemTypePath $dir
	
	
	
	#------------------
	set file [open c:/tcl/lib/pepe.txt w 0600] 
    puts $file "esto es una prueba" 
    close $file
	#------------------
}

#-------------------------------------------------------------------------------


proc AfterReadGIDProject { filename } {
    
	
	#------------------
	set file [open c:/tcl/lib/pepe.txt w 0600] 
    puts $file "esto es una prueba2" 
    close $file
	#------------------
	
	
	
    # Save ProblemPath
    set projectpath $filename
    append projectpath .gid
    set ::FemDemKratos::ProblemPath $projectpath
    
	
#	#------------------
#	set file [open c:/tcl/lib/pepe.txt w 0600] 
#    puts $file "------------" 
#	puts $file $projectpath 
#    close $file
#	#------------------
	
	
    # Save ProblemName
    # if {$::tcl_platform(platform) eq "windows"} {}
    if {[regexp -all {\\} $filename] > 0} {
        # Windows
        regsub -all {\\} $filename { } filename
    } else {
        # Unix
        regsub -all {/} $filename { } filename
    }
    set filename [lreplace $filename 0 [expr { [llength $filename]-2 }]]
    set ::FemDemKratos::ProblemName $filename
}



proc BeforeRunCalculation { batfilename basename dir problemtypedir gidexe args } {  
    # Set Parallel Configuration
    set paralleltype [GiD_AccessValue get gendata Parallel_Configuration]
    
    # Write Initial fractures data
    #if {([GiD_AccessValue get gendata Fracture_Propagation] eq true) && ($paralleltype ne "MPI")} {
    #    # Define GiDPath
    #    if {[regexp -all {\\} $gidexe] > 0} {
    #        # Windows
    #        regsub -all {\\} $gidexe {/} gidexe
    #    }
    #    set gidexe [string trimright $gidexe gid.exe]
    #    
    #    if {[GiD_AccessValue get gendata Domain_Size] eq 2} {
    #    
    #        source [file join $problemtypedir FracturePropagation2D.tcl]
    #        WriteInitialFracturesData $dir $problemtypedir $gidexe
    #    
    #    } else {
    #        
    #        source [file join $problemtypedir FracturePropagation3D.tcl]
    #        WriteInitialFracturesData $dir $problemtypedir $gidexe
    #    }
    #}

    # Write MDPA
    source [file join $problemtypedir Mdpa.tcl]
    set TableDict [WriteMdpa $basename $dir $problemtypedir]

    # Write ProjectParameters
    source [file join $problemtypedir ProjectParameters.tcl]
    WriteProjectParameters $basename $dir $problemtypedir $TableDict
    
     # # Copy python script in the problemdir
     # if {[GiD_AccessValue get gendata Fracture_Propagation] eq true} {
     #     file copy -force [file join $problemtypedir poromechanics_fracture_main.py] [file join $dir MainKratos.py]
     # } else {
     #     file copy -force [file join $problemtypedir poromechanics_main.py] [file join $dir MainKratos.py]
     # }
    
    # Run the problem
      #set run 1
      #catch {
      #    if {$paralleltype eq "MPI"} {set run 0}
      #}
      #if {$run} {
      #    return ""
      #} {
      #    return [list "-cancel-" [= "You have selected MPI parallelism system.\n\
      #                                Input files have been written.\n\
      #                                Run the case with: mpirun -np \[npartitions\] python3 MainKratos.py" ]]
      #}
    
    ### Measure time
    #set start_time_1 [clock clicks]
    #set end_time_1 [expr { [clock clicks]-$start_time_1 }]
    #WarnWin "Time for GenerateNewFractures: $end_time_1 clicks"
    ###
}