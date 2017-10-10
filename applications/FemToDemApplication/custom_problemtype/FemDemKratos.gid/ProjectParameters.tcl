proc WriteProjectParameters { basename dir problemtypedir TableDict} {

    ## Source auxiliar procedures
    source [file join $problemtypedir ProjectParametersAuxProcs.tcl]
        
    ## Start ProjectParameters.json file
    set filename [file join $dir ProjectParameters.json]
    set FileVar [open $filename w]
    
    puts $FileVar "\{"

    ## AMR data
    puts $FileVar "   \"AMR_data\": \{"
    puts $FileVar "        \"activate_AMR\":                    [GiD_AccessValue get gendata Activate_AMR],"
    puts $FileVar "        \"plane_state\":                    \"[GiD_AccessValue get gendata Plane_state]\","
    puts $FileVar "        \"mesh_optimality_criteria\":       \"[GiD_AccessValue get gendata Mesh_Optimality_Criteria]\","
    puts $FileVar "        \"permissible_error\":               [GiD_AccessValue get gendata Permissible_Error],"
    puts $FileVar "        \"refinement_frequency\":            [GiD_AccessValue get gendata Refinement_Frequency],"
    puts $FileVar "        \"gid_path\":                       \"[GiD_AccessValue get gendata gid_path]\","
    puts $FileVar "    \},"
    ## problem_data
    puts $FileVar "   \"problem_data\": \{"
    puts $FileVar "        \"problem_name\":         \"$basename\","
    puts $FileVar "        \"model_part_name\":      \"Structure\","
    puts $FileVar "        \"domain_size\":          [GiD_AccessValue get gendata Domain_Size],"
    puts $FileVar "        \"start_time\":           [GiD_AccessValue get gendata Start_Time],"
    puts $FileVar "        \"end_time\":             [GiD_AccessValue get gendata End_Time],"
    puts $FileVar "        \"time_step\":            [GiD_AccessValue get gendata Delta_Time],"
	puts $FileVar "        \"echo_level\":           [GiD_AccessValue get gendata Echo_Level],"
    puts $FileVar "    \},"
    ## solver_settings
    puts $FileVar "   \"solver_settings\": \{"
    if {[GiD_AccessValue get gendata Parallel_Configuration] eq "MPI"} {
        puts $FileVar "        \"solver_type\":                        \"solid_mechanics_implicit_dynamic_solver\","
    } else {
        puts $FileVar "        \"solver_type\":                        \"solid_mechanics_implicit_dynamic_solver\","
    }
	puts $FileVar "            \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
	puts $FileVar "            \"solution_type\":                      \"Dynamic\","
	puts $FileVar "            \"time_integration_method\":            \"Newmark\","	
	

    puts $FileVar "            \"model_import_settings\":              \{"
    puts $FileVar "                 \"input_type\":         \"mdpa\","
    puts $FileVar "                 \"input_filename\":     \"$basename\","
    puts $FileVar "                 \"input_file_label\":    0"
    puts $FileVar "            \},"
   # puts $FileVar "        \"buffer_size\":                        2,"
   # puts $FileVar "        \"echo_level\":                         [GiD_AccessValue get gendata Echo_Level],"
   # puts $FileVar "        \"clear_storage\":                      false,"
   # puts $FileVar "        \"compute_reactions\":                  [GiD_AccessValue get gendata Write_Reactions],"
   # puts $FileVar "        \"move_mesh_flag\":                     [GiD_AccessValue get gendata Move_Mesh],"
   # set IsPeriodic [GiD_AccessValue get gendata Periodic_Interface_Conditions]
   puts $FileVar "             \"line_search\":                          false,"
   puts $FileVar "             \"convergence_criterion\":               \"[GiD_AccessValue get gendata Convergence_Criterion]\","
   puts $FileVar "             \"displacement_relative_tolerance\":      [GiD_AccessValue get gendata Displacement_Relative_Tolerance],"
   puts $FileVar "             \"displacement_absolute_tolerance\":      [GiD_AccessValue get gendata Displacement_Absolute_Tolerance],"
   puts $FileVar "             \"residual_relative_tolerance\":          [GiD_AccessValue get gendata Residual_Relative_Tolerance],"
   puts $FileVar "             \"residual_absolute_tolerance\":          [GiD_AccessValue get gendata Residual_Absolute_Tolerance],"
   puts $FileVar "             \"max_iteration\":                        [GiD_AccessValue get gendata Max_Iterations],"

    puts $FileVar "        \"linear_solver_settings\":     \{"
            puts $FileVar "          \"solver_type\":   \"SuperLUSolver\""
			puts $FileVar "          \"scaling\":       \"false\""

    puts $FileVar "        \},"

	set PutStrings \[

    set BGroups [GiD_Info conditions Body_Part groups]
    #W "tamanyo [llength $Groups]"

    # Body_Part
    if {[llength $BGroups] eq "1"} {
    AppendGroupName PutStrings Body_Part
    } else {
    AppendGroupNames PutStrings Body_Part
    }
    
    append PutStrings \]
    puts $FileVar "        \"problem_domain_sub_model_part_list\": $PutStrings,"
    ## processes_sub_model_part_list
    set PutStrings \[
    # Solid_Displacement
    AppendGroupNames PutStrings Solid_Displacement
    # Force
    AppendGroupNames PutStrings Force
    # Face_Load
    AppendGroupNames PutStrings Face_Load
    # Normal_Load
    AppendGroupNames PutStrings Normal_Load
    # Body_Acceleration
    AppendGroupNames PutStrings Body_Acceleration
    set PutStrings [string trimright $PutStrings ,]
    append PutStrings \]
    puts $FileVar "        \"processes_sub_model_part_list\":      $PutStrings,"
    ## body_domain_sub_model_part_list
    set PutStrings \[
    AppendGroupNames PutStrings Body_Part


    puts $FileVar "     \},"
    











    puts $FileVar ""
    puts $FileVar "\}"

    close $FileVar
}
