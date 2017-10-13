from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FemToDemApplication  import *
import shutil 
CheckForPreviousImport()

def Wait():
    input("Press Something")

#============================================================================================================================
class AdaptiveMeshRefinementUtility:

    def __init__(self, ProjectParameters, starting_time, solver_constructor, constitutive_law_utility,
     gid_output_utility, conditions_util, ProblemPath):
        
        ## Parameters and utilities
        self.ProjectParameters = ProjectParameters
        self.solver_constructor = solver_constructor
        self.constitutive_law_utility = constitutive_law_utility
        self.gid_output_utility = gid_output_utility
        
        self.conditions_util = conditions_util
        self.problem_path    = os.getcwd()
        self.AMR_files_path  = os.path.join(self.problem_path, "AMR_Files")
        


        #print("path", self.problem_path)
        #print("path amr!", self.AMR_files_path)
        #Wait()
        self.n_refinements = 0
        self.last_refinement_id = 1 #Same as the initial current_id
        
        ## Time operations initialization
        self.ending_time = ProjectParameters["problem_data"]["end_time" ].GetDouble()
        self.delta_time  = ProjectParameters["problem_data"]["time_step"].GetDouble()
        
        # set AMR frequency
        self.amr_frequency = ProjectParameters["AMR_data"]["refinement_frequency"].GetDouble()
        if(self.amr_frequency < self.delta_time):
            self.amr_frequency = self.delta_time

        # set time counter
        self.time_counter = starting_time + self.amr_frequency

        # set time operation tolerance
        self.tolerance = self.delta_time * 1e-10;

#============================================================================================================================
    def Initialize(self):
        
        self.gid_path = self.ProjectParameters["AMR_data"]["gid_path"].GetString()
        AMR_info_path = os.path.join(self.problem_path,"AMR_info.txt")

        #print("gid", self.gid_path)
        #print(len(str(self.gid_path)))
        #print("lera",str(self.gid_path)[0])
        #self.gid_path = str(self.gid_path)[:-10]
        #print("gid", self.gid_path)
        Wait()
        activate_AMR = True
        
        if(self.amr_frequency > self.ending_time):
            activate_AMR = False
        
        # creates the AMR_Files folder or roemove if exists
        if not os.path.isdir(self.AMR_files_path):
            os.makedirs(str(self.AMR_files_path))
        else:  
            shutil.rmtree(str(self.AMR_files_path), ignore_errors = True)
            os.makedirs(str(self.AMR_files_path))

        if os.path.isfile(str(AMR_info_path)):
            shutil.rmtree(str(AMR_info_path), ignore_errors = True)
        
        return activate_AMR
        
#============================================================================================================================
    def CheckAMR(self, current_time):
        
        refine = False
        last_mesh = False
        
        if( (current_time + self.tolerance >= self.time_counter) and (current_time + self.tolerance < self.ending_time) ):
            self.time_counter = self.time_counter + self.amr_frequency
            refine = True
        elif(current_time + self.tolerance >= self.ending_time):
            last_mesh = True

        return refine, last_mesh
        
#============================================================================================================================       
    def Execute(self, model_part, main_step_solver, gid_output_util, current_time, current_id):
        
        ## Previous definitions ---------------------------------------------------------------------------------------------
        
        problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()
        #GidOutputConfiguration = self.ProjectParameters.GidOutputConfiguration
        output_mode = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()
        output_multiple_files = self.ProjectParameters["output_configuration"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()
        plane_state = self.ProjectParameters["AMR_data"]["plane_state"].GetString()
        mesh_optimality_criteria = self.ProjectParameters["AMR_data"]["mesh_optimality_criteria"].GetString()
        permissible_error = self.ProjectParameters["AMR_data"]["permissible_error"].GetDouble()

        ## Finalize previous post results -----------------------------------------------------------------------------------
        #print("dentro de execute1", os.getcwd())
        #i = 69
        #print(str("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "))
        #Wait()
        gid_output_util.finalize_results()
        
        if(output_mode=="GiD_PostBinary"):
            if(output_multiple_files=="MultipleFiles"):
                for i in range(self.last_refinement_id,current_id + 1):
                    # os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.bin")
                    src = os.path.join(self.problem_path, str(problem_name) + str(i) + ".post.bin")
                    dst = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) +"_step_" + str(i) + ".post.bin")
                    shutil.move(src, dst)
            else:
                # os.system("move "+str(self.problem_path)+"/"+str(problem_name)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.bin")
                src = os.path.join(self.problem_path, str(problem_name) + ".post.bin ")
                dst = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.bin")
                shutil.move(src, dst)

        else:
            if(output_multiple_files=="MultipleFiles"):
                for i in range(self.last_refinement_id,current_id+1):
                    #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.msh")
                    #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.res")
                    src_mesh = os.path.join(self.problem_path, str(problem_name) + "_" + str(i) + ".post.msh")
                    src_res  = os.path.join(self.problem_path, str(problem_name) + "_" + str(i) + ".post.res")
                    dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_"+  str(self.n_refinements) + "_step_" + str(i) + ".post.msh")
                    dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + "_step_" + str(i) + ".post.res")

                    shutil.move(src_mesh, dst_mesh)
                    shutil.move(src_res , dst_res)

            else:
                #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.msh")
                #os.system("move "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.res")

                src_mesh = os.path.join(self.problem_path, str(problem_name) + "_" + str(self.last_refinement_id) +".post.msh")
                src_res  = os.path.join(self.problem_path, str(problem_name) + "_" + str(self.last_refinement_id) +".post.res" )
                dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.msh")
                dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_results_mesh_" + str(self.n_refinements) + ".post.res")

                shutil.move(src_mesh, dst_mesh)
                shutil.move(src_res , dst_res)

        #print("antes de copy  ")   #cornejo
        #Wait()        
        #os.system("copy "+str(self.problem_path)+"/"+str(problem_name)+".mdpa "+str(self.AMR_files_path)+"/"+str(problem_name)+"_mesh_"+str(self.n_refinements)+".mdpa")
        src_mdpa = os.path.join(self.problem_path, str(problem_name) +".mdpa")
        dst_mdpa = os.path.join(self.AMR_files_path, str(problem_name) + "_mesh_" + str(self.n_refinements) + ".mdpa")
        shutil.copy(src_mdpa, dst_mdpa)


        ## MESH LOOP --------------------------------------------------------------------------------------------------------
        
        #print("dentro de execute2")   #cornejo
        #Wait()
        
        print("----------------------------------")
        print(" START MESH REFINEMENT ITERATIONS ")
        print("----------------------------------")
        
        mesh_convergence = False
        max_num_iter = 5
        iteration_number = 0
        
        while(mesh_convergence == False and iteration_number < max_num_iter ):
            
            iteration_number = iteration_number + 1
            print("MESH ITERATION : ", iteration_number)
            
            ## Generate files for new mesh ----------------------------------------------------------------------------------
            #print("aqui peta")
            #Wait()

            if(iteration_number == 1):
                AdaptiveMeshRefinementProcess(model_part,plane_state,
                                                         problem_name,
                                                         self.problem_path,
                                                         mesh_optimality_criteria,
                                                         permissible_error,
                                                         self.n_refinements).Execute()
                Wait()

                # Move the posts of the amr to the AMR_folder
                src_mesh = os.path.join(self.problem_path, str(problem_name)   + "_AMR_parameters.post.msh")
                src_res  = os.path.join(self.problem_path, str(problem_name)   + "_AMR_parameters.post.res")
                dst_mesh = os.path.join(self.AMR_files_path, str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.msh")
                dst_res  = os.path.join(self.AMR_files_path, str(problem_name) + "_AMR_parameters_mesh_" + str(self.n_refinements) + ".post.res")

                shutil.move(src_mesh, dst_mesh)
                shutil.move(src_res, dst_res)

            else:
                AdaptiveMeshRefinementProcess(model_part, plane_state,
                                                          problem_name,
                                                          self.problem_path,
                                                          mesh_optimality_criteria,
                                                          permissible_error,
                                                          self.n_refinements).ExecuteAfterOutputStep()
            
            print("despues execute after output step")
            Wait()
            

            # GID GENERATE THE NEW MESH BASED ON THE BACKGROUND MESH  ----> TODO
            #os.system("cd " + str(self.gid_path) + " && ./gid -b "+ str(self.problem_path)+"/"+str(problem_name)+".bch -n")
            #os.system("cd && mv "+str(self.problem_path)+"/"+str(problem_name)+".dat "+str(self.problem_path)+"/"+str(problem_name)+".mdpa && rm "+str(self.problem_path)+"/"+str(problem_name)+"-1.dat")
            
            # Execute .bch file with GiD 
            os.system("cd " + str(self.gid_path[:-8]) + " && gid -b " + os.path.join(self.problem_path, str(problem_name) + ".bch -n"))
            print("despues de gid commands gid path: ", self.gid_path[:-8])
            Wait()            
            #shutil.move(os.path.join(str(self.problem_path), str(problem_name) + ".dat" ) , os.path.join(str(self.problem_path), str(problem_name) + ".mdpa"))
            #shutil.rmtree(os.path.join(str(self.problem_path), str(problem_name)+"-1.dat"))
            
            print("despues de removees")
            Wait()
            #------------------------------------------------------------------->> Aqui estamos
            
            ## Finalize previous mesh ---------------------------------------------------------------------------------------
            
            #Finalize previous solver
            main_step_solver.Finalize()
            
            #Save previous model_part
            model_part_old = model_part
            
            ## Generate new Model Part --------------------------------------------------------------------------------------

            # Definition of model part
            model_part = ModelPart("SolidDomain")

            # Set problem variables
            self.solver_constructor.AddVariables(model_part)

            # Reading model part
            model_part_io = ModelPartIO(problem_name)
            model_part_io.ReadModelPart(model_part)

            # Set buffer size
            buffer_size = 2
            model_part.SetBufferSize(buffer_size)

            # Set degrees of freedom
            self.solver_constructor.AddDofs(model_part)
            
            # Set ProcessInfo variables and fill the previous steps of the buffer with the initial conditions
            current_time = current_time - (buffer_size-1)*self.delta_time
            model_part.ProcessInfo[TIME] = current_time
            for step in range(buffer_size-1):
                current_time = current_time + self.delta_time
                model_part.CloneTimeStep(current_time)

            ## Initialize new model part ------------------------------------------------------------------------------------
            
            # Definition of output utility
            gid_output_util = self.gid_output_utility.GidOutputUtility(GidOutputConfiguration,problem_name,current_time,self.ending_time,self.delta_time)
            
            # Set constitutive laws
            self.constitutive_law_utility.SetConstitutiveLaw(model_part)

            # Define and initialize the main solver
            main_step_solver = self.solver_constructor.CreateSolver(model_part, self.ProjectParameters.SolverSettings)
            model_part.ProcessInfo[MESH_REFINED] = 1
            main_step_solver.Initialize()

            # Initialize imposed conditions
            NonLinearImposedConditions = self.conditions_util.Initialize(model_part)
            
            ## Mapping of variables -----------------------------------------------------------------------------------------
            
            MappingVariablesProcess(model_part_old,model_part,self.ProjectParameters.ConditionsOptions.Imposed_Displacement).Execute()

            ## Erase old Model Part -----------------------------------------------------------------------------------------
            
            model_part_old = None
            
            ## Test new Mesh ------------------------------------------------------------------------------------------------
            
            mesh_convergence = main_step_solver.TestNewMesh()
            
        if(mesh_convergence==True):
            print("NEW MESH CONVERGED AFTER ",iteration_number," ITERATIONS")
        else:
            print("### WARNING: NO MESH CONVERGED AFTER ", iteration_number, " ITERATIONS ###")
        
        print(model_part)
        
        ## Saving files of new mesh -----------------------------------------------------------------------------------------
        
        os.system("cp "+str(self.problem_path)+"/"+str(problem_name)+".bgm "+str(self.AMR_files_path)+"/"+str(problem_name)+"_AMR_"+str(self.n_refinements+1)+".bgm")
        
        ## Initialize post files of new mesh --------------------------------------------------------------------------------

        gid_output_util.initialize_results(model_part, current_id+1) #For single post file
        
        self.last_refinement_id = current_id+1
        self.n_refinements = self.n_refinements + 1
        
        return model_part, main_step_solver, gid_output_util

#============================================================================================================================
    def Finalize(self, model_part, current_id):
        
        ## Previous definitions ---------------------------------------------------------------------------------------------
        
        problem_name = self.ProjectParameters.problem_name
        GidOutputConfiguration = self.ProjectParameters.GidOutputConfiguration
        output_mode = GidOutputConfiguration.GiDPostMode
        output_multiple_files = GidOutputConfiguration.GiDPostFiles
        plane_state = self.ProjectParameters.plane_state
        mesh_optimality_criteria = self.ProjectParameters.mesh_optimality_criteria
        permissible_error = self.ProjectParameters.permissible_error

        ## Finalize previous post results -----------------------------------------------------------------------------------
        
        if(output_mode=="Binary"):
            if(output_multiple_files=="Multiples"):
                for i in range(self.last_refinement_id,current_id+1):
                    os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.bin")
            else:
                os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+".post.bin "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.bin")
        else:
            if(output_multiple_files=="Multiples"):
                for i in range(self.last_refinement_id,current_id+1):
                    os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.msh")
                    os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(i)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+"_step_"+str(i)+".post.res")
            else:
                os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.msh")
                os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_"+str(self.last_refinement_id)+".post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_results_mesh_"+str(self.n_refinements)+".post.res")
        
        os.system("cp "+str(self.problem_path)+"/"+str(problem_name)+".mdpa "+str(self.AMR_files_path)+"/"+str(problem_name)+"_mesh_"+str(self.n_refinements)+".mdpa")
        
        ## Compute and save info of last mesh -------------------------------------------------------------------------------
        
        AdaptiveMeshRefinementProcess(model_part,plane_state,problem_name,self.problem_path,mesh_optimality_criteria,permissible_error,self.n_refinements).ExecuteFinalize()
        
        os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_AMR_parameters.post.msh "+str(self.AMR_files_path)+"/"+str(problem_name)+"_AMR_parameters_mesh_"+str(self.n_refinements)+".post.msh")
        os.system("mv "+str(self.problem_path)+"/"+str(problem_name)+"_AMR_parameters.post.res "+str(self.AMR_files_path)+"/"+str(problem_name)+"_AMR_parameters_mesh_"+str(self.n_refinements)+".post.res")

        os.system("mv "+str(self.problem_path)+"/AMR_info.txt "+str(self.AMR_files_path))
#============================================================================================================================        
