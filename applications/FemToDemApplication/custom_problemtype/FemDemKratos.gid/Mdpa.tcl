proc WriteMdpa { basename dir problemtypedir } {
    
    ## Source auxiliar procedures
    source [file join $problemtypedir MdpaAuxProcs.tcl]
    
    ## Start MDPA file
    set filename [file join $dir ${basename}.mdpa]
    set FileVar [open $filename w]
    
	set TableId 0
    set TableDict [dict create]
    # Solid_Displacement
    ConstraintVectorTable FileVar TableId TableDict Solid_Displacement DISPLACEMENT
    # Fluid_Pressure
    PressureTable FileVar TableId TableDict Fluid_Pressure WATER_PRESSURE
    # Force
    VectorTable FileVar TableId TableDict Force FORCE
    # Face_Load
    VectorTable FileVar TableId TableDict Face_Load FACE_LOAD
    # Normal_Load
    # Body_Acceleration
    VectorTable FileVar TableId TableDict Body_Acceleration VOLUME_ACCELERATION
    puts $FileVar ""
    
    ## Properties
    set PropertyId 0
    set PropertyDict [dict create]
	
	set Groups [GiD_Info conditions Body_Part groups]
	for {set i 0} {$i < [llength $Groups]} {incr i} {
		
		puts $FileVar "Begin Properties $PropertyId"
		puts $FileVar "  YIELD_SURFACE [lindex [lindex $Groups $i] 3]"
		puts $FileVar "  YOUNG_MODULUS [lindex [lindex $Groups $i] 4]"
		puts $FileVar "  DENSITY [lindex [lindex $Groups $i] 5]"
		puts $FileVar "  POISSON_RATIO [lindex [lindex $Groups $i] 6]"
		puts $FileVar "  THICKNESS [lindex [lindex $Groups $i] 7]"
		puts $FileVar "  YIELD_STRESS_C [lindex [lindex $Groups $i] 8]"
		puts $FileVar "  YIELD_STRESS_T [lindex [lindex $Groups $i] 9]"
		puts $FileVar "  FRAC_ENERGY_T [lindex [lindex $Groups $i] 10]"
		puts $FileVar "  INTERNAL_FRICTION_ANGLE [lindex [lindex $Groups $i] 11]"
		puts $FileVar "  RAYLEIGH_BETA [lindex [lindex $Groups $i] 12]"
		puts $FileVar "  RAYLEIGH_ALPHA [lindex [lindex $Groups $i] 13]"
		puts $FileVar "End Properties"
		puts $FileVar ""
	}
	
	## Nodes
    set Nodes [GiD_Info Mesh Nodes]
    puts $FileVar "Begin Nodes"
    for {set i 0} {$i < [llength $Nodes]} {incr i 4} {
        puts $FileVar "  [lindex $Nodes $i]  [lindex $Nodes [expr { $i+1 }]] [lindex $Nodes [expr { $i+2 }]] [lindex $Nodes [expr { $i+3 }]]"
        #puts -nonewline $FileVar "  [lindex $Nodes $i]  "
        #puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+1 }]]]
        #puts -nonewline $FileVar " "
        #puts -nonewline $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+2 }]]]
        #puts -nonewline $FileVar " "
        #puts $FileVar [format  "%.10f" [lindex $Nodes [expr { $i+3 }]]]
    }
    puts $FileVar "End Nodes"
    puts $FileVar ""
    puts $FileVar ""
    
	## Elements
	set FIC [GiD_AccessValue get gendata FIC_Stabilization]
    set IsQuadratic [GiD_Info Project Quadratic]
	# UPwSmallStrainElement2D3N
	for {set i 0} {$i < [llength $Groups]} {incr i} {
		WriteElements FileVar [lindex $Groups $i] triangle AleCornVelElement $BodyElemsProp Triangle2D3Connectivities
	}
	
	## TODO Conditions line 365 problem de POROMECH
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}