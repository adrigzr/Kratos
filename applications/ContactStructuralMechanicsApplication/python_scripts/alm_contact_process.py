from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ALMContactProcess(Model, settings["Parameters"])

import python_process
##all the processes python processes should be derived from "python_process"
class ALMContactProcess(python_process.PythonProcess):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "contact_model_part"          : "Contact_Part",
            "axisymmetric"                : false,
            "assume_master_slave"         : "",
            "contact_type"                : "Frictionless",
            "frictional_law"              : "Coulomb",
            "search_factor"               : 2.0,
            "active_check_factor"         : 0.01,
            "dual_search_check"           : false,
            "strict_search_check"         : true,
            "max_number_results"          : 1000,
            "bucket_size"                 : 4,
            "normal_variation"            : false,
            "pair_variation"              : true,
            "manual_ALM"                  : false,
            "stiffness_factor"            : 10.0,
            "penalty_scale_factor"        : 1.0,
            "use_scale_factor"            : true,
            "penalty"                     : 0.0,
            "scale_factor"                : 1.0e0,
            "tangent_factor"              : 0.1,
            "type_search"                 : "InRadius",
            "use_exact_integration"       : true,
            "hard_clear_after_step"       : false,
            "database_step_update"        : 1,
            "integration_order"           : 2,
            "predict_with_linear_solver"  : false,
            "max_gap_factor"              : 0.0,
            "linear_solver_settings"      : {
                "solver_type"             : "SuperLUSolver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 1
            },
            "debug_mode"                  : false
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
   
        self.main_model_part = model_part[self.params["model_part_name"].GetString()]
        self.computing_model_part_name = self.params["computing_model_part_name"].GetString()

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.contact_model_part = model_part[self.params["contact_model_part"].GetString()]
        
        if (self.params["assume_master_slave"].GetString() != ""):
            for node in self.contact_model_part.Nodes:
                node.Set(KratosMultiphysics.SLAVE, False)
            del(node)
            model_part_slave = self.main_model_part.GetSubModelPart(self.params["assume_master_slave"].GetString())
            for node in model_part_slave.Nodes:
                node.Set(KratosMultiphysics.SLAVE, True)
            del(node)
        
        self.axisymmetric  = self.params["axisymmetric"].GetBool()
        if (self.axisymmetric == True) and (self.dimension == 3):
            raise NameError("3D and axisymmetric makes no sense")
        self.normal_variation = self.params["normal_variation"].GetBool()
        self.database_step_update = self.params["database_step_update"].GetInt() 
        self.database_step = 0
        self.frictional_law = self.params["frictional_law"].GetString()
        self.predict_with_linear_solver = self.params["predict_with_linear_solver"].GetBool()
        self.debug_mode = self.params["debug_mode"].GetBool()
        self.hard_clear_after_step = self.params["hard_clear_after_step"].GetBool()
        
        # Debug
        if (self.debug_mode == True):
            self.output_file = "POSTSEARCH"

            self.gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            self.singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
            self.deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            self.write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
        
    def ExecuteInitialize(self):
        # Appending the conditions created to the self.main_model_part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")

        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.params["contact_type"].GetString() == "Frictional":
            computing_model_part.Set(KratosMultiphysics.SLIP, True) 
        else:
            computing_model_part.Set(KratosMultiphysics.SLIP, False) 
            
        # We recompute the normal at each iteration (false by default)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.CONSIDER_NORMAL_VARIATION] = self.normal_variation
        # We recompute the pairs at each iteration (true by default)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.CONSIDER_PAIR_VARIATION] = self.params["pair_variation"].GetBool()
        # We set the max gap factor for the gap adaptation
        max_gap_factor = self.params["max_gap_factor"].GetDouble()
        if (max_gap_factor > 0.0):
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ADAPT_PENALTY] = True
        else:
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ADAPT_PENALTY] = False
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.MAX_GAP_FACTOR] = max_gap_factor
        
        # We set the value that scales in the tangent direction the penalty and scale parameter
        if self.params["contact_type"].GetString() == "Frictional":
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.TANGENT_FACTOR] = self.params["tangent_factor"].GetDouble()
        
        # Copying the properties in the contact model part
        self.contact_model_part.SetProperties(computing_model_part.GetProperties())
        
        # Setting the integration order and active check factor
        for prop in computing_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.params["integration_order"].GetInt() 
            prop[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = self.params["active_check_factor"].GetDouble()
            
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(computing_model_part)
        
        if self.params["contact_type"].GetString() == "Frictionless":
            if self.normal_variation == True:
                if self.axisymmetric == True:
                    condition_name = "ALMNVFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionlessMortarContact"
            else:
                if self.axisymmetric == True:
                    condition_name = "ALMFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionlessMortarContact"
        elif self.params["contact_type"].GetString() == "Frictional":
            if self.normal_variation == True:
                if self.axisymmetric == True:
                    condition_name = "ALMNVFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionalMortarContact"
            else:
                if self.axisymmetric == True:
                    condition_name = "ALMFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionalMortarContact"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(computing_model_part) 
        
        # It should create the conditions automatically
        interface_parameters = KratosMultiphysics.Parameters("""{"condition_name": "", "final_string": "", "simplify_geometry": false}""")
        interface_parameters["condition_name"].SetString(condition_name)
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(computing_model_part, self.contact_model_part, interface_parameters) 
        else:
            self.Preprocess.GenerateInterfacePart3D(computing_model_part, self.contact_model_part, interface_parameters) 

        # When all conditions are simultaneously master and slave
        if (self.params["assume_master_slave"].GetString() == ""):
            for node in self.contact_model_part.Conditions:
                node.Set(KratosMultiphysics.SLAVE, True)
            del(node)
            for cond in self.contact_model_part.Conditions:
                cond.Set(KratosMultiphysics.SLAVE, True)
            del(cond)

        if (self.params["manual_ALM"].GetBool() == False):
            # Computing the scale factors or the penalty parameters (StiffnessFactor * E_mean/h_mean)
            self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(computing_model_part)
            self.find_nodal_h.Execute()
            
            alm_var_parameters = KratosMultiphysics.Parameters("""{}""")
            alm_var_parameters.AddValue("stiffness_factor",self.params["stiffness_factor"])
            alm_var_parameters.AddValue("penalty_scale_factor",self.params["penalty_scale_factor"])
            self.alm_var_process = ContactStructuralMechanicsApplication.ALMVariablesCalculationProcess(self.contact_model_part,KratosMultiphysics.NODAL_H, alm_var_parameters)
            self.alm_var_process.Execute()
            # We don't consider scale factor
            if (self.params["use_scale_factor"].GetBool() == False):
                self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR] = 1.0
        else:
            # We set the values in the process info
            self.main_model_part.ProcessInfo[KratosMultiphysics.INITIAL_PENALTY] = self.params["penalty"].GetDouble()
            self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR] = self.params["scale_factor"].GetDouble()
            
        # We print the parameters considered
        print("The parameters considered finally are: ")            
        print("SCALE_FACTOR: ","{:.2e}".format(self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR]))
        print("INITIAL_PENALTY: ","{:.2e}".format(self.main_model_part.ProcessInfo[KratosMultiphysics.INITIAL_PENALTY]))
            
        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(computing_model_part)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)

        # Creating the search
        search_parameters = KratosMultiphysics.Parameters("""{}""")
        search_parameters.AddValue("type_search",self.params["type_search"])
        search_parameters.AddValue("allocation_size",self.params["max_number_results"])
        search_parameters.AddValue("bucket_size",self.params["bucket_size"])
        search_parameters.AddValue("search_factor",self.params["search_factor"])
        search_parameters.AddValue("dual_search_check",self.params["dual_search_check"])
        search_parameters.AddValue("strict_search_check",self.params["strict_search_check"])
        search_parameters.AddValue("use_exact_integration",self.params["use_exact_integration"])
        self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch(computing_model_part, search_parameters)
        
        # We initialize the conditions    
        self.alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.contact_model_part) 
        self.alm_init_var.Execute()
        
        # We initialize the search utility
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()
        
    def ExecuteBeforeSolutionLoop(self):
        if self.params["contact_type"].GetString() == "Frictionless":  
            self.contact_search.TotalClearALMFrictionlessMortarConditions()
        else:
            self.contact_search.TotalClearComponentsMortarConditions()
    
    def ExecuteInitializeSolutionStep(self):
        self.database_step += 1
        self.global_step = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME_STEPS]
        
        if (self.database_step >= self.database_step_update or self.global_step == 1):
            # We solve one linear step with a linear strategy if needed
            
            self.contact_search.UpdateMortarConditions()
            #self.contact_search.CheckMortarConditions()
                
            if (self.predict_with_linear_solver == True and self.global_step > 1):
                # Debug
                if (self.debug_mode == True):
                    self._debug_output(self.global_step, "_LINEARPRED")
                self._linear_solver_predict()
                if (self.hard_clear_after_step == True):
                    if self.params["contact_type"].GetString() == "Frictionless":  
                        self.contact_search.TotalClearALMFrictionlessMortarConditions()
                    else:
                        self.contact_search.TotalClearComponentsMortarConditions()
                        
                    self.contact_search.UpdateMortarConditions()
                    #self.contact_search.CheckMortarConditions()
                else:
                    self.contact_search.CleanMortarConditions()
                
            # Debug
            if (self.debug_mode == True):
               self._debug_output(self.global_step, "")
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if (self.database_step >= self.database_step_update or self.global_step == 1):
            self._clear_sets()
            
    def ExecuteFinalize(self):
        pass

    def _clear_sets(self):
        if (self.hard_clear_after_step == True):
            if self.params["contact_type"].GetString() == "Frictionless":  
                self.contact_search.TotalClearALMFrictionlessMortarConditions()
            else:
                self.contact_search.TotalClearComponentsMortarConditions()
        else:
            if self.params["contact_type"].GetString() == "Frictionless":
                self.contact_search.PartialClearALMFrictionlessMortarConditions()
            else:
                self.contact_search.PartialClearComponentsMortarConditions()
                
        self.database_step = 0
                
    def _linear_solver_predict(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.params["linear_solver_settings"])
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part, 
                                                                scheme, 
                                                                linear_solver, 
                                                                builder_and_solver, 
                                                                compute_reactions, 
                                                                reform_step_dofs, 
                                                                calculate_norm_dx,
                                                                move_mesh_flag
                                                                )
        strategy.SetEchoLevel(0)
        strategy.Check()
        strategy.Solve()
        strategy.Clear()
    
    def _debug_output(self, label, name):

        gid_io = KratosMultiphysics.GidIO(self.output_file+name+"_STEP_"+str(label), self.gid_mode, self.singlefile, self.deformed_mesh_flag, self.write_conditions)
        
        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.main_model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.main_model_part.GetMesh())
        
        gid_io.WriteNodalFlags(KratosMultiphysics.INTERFACE, "INTERFACE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KratosMultiphysics.ACTIVE, "ACTIVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KratosMultiphysics.SLAVE, "SLAVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(KratosMultiphysics.NORMAL, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResults(KratosMultiphysics.NORMAL_CONTACT_STRESS, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResults(ContactStructuralMechanicsApplication.WEIGHTED_GAP, self.main_model_part.Nodes, label, 0)
        
        gid_io.FinalizeResults()
        
        #raise NameError("DEBUG")
