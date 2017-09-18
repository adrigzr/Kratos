from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import system python modules
import time as timer
import os

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication     as KratosSolid
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.ExternalSolversApplication    as KratosSolvers
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import MainSolidFEM
import main_script as MainDEM

def Wait():
	input("Press Something")


class FEM_Solution(MainSolidFEM.Solution):

	def Info(self):
		print("FEM part of the FEMDEM application")

	                
	def __init__(self):

		#### TIME MONITORING START ####
		# Time control starts        
		print(timer.ctime())
		# Measure process time
		self.t0p = timer.clock()
		# Measure wall time
		self.t0w = timer.time()
		#### TIME MONITORING END ####


		#### PARSING THE PARAMETERS ####

		# Import input
		parameter_file = open("ProjectParameters.json",'r')
		self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

		# set echo level
		self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()

		print(" ")

		# defining the number of threads:
		num_threads =  self.GetParallelSize()
		print("::[KSM Simulation]:: [OMP USING",num_threads,"THREADS ]")
		#parallel.PrintOMPInfo()


		print(" ")
		print("::[KSM Simulation]:: [Time Step:", self.ProjectParameters["problem_data"]["time_step"].GetDouble()," echo:", self.echo_level,"]")

		#### Model_part settings start ####

		# Defining the model_part
		self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

		#self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DIMENSION, self.ProjectParameters["problem_data"]["dimension"].GetInt())
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, self.ProjectParameters["problem_data"]["time_step"].GetDouble())
		self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, self.ProjectParameters["problem_data"]["start_time"].GetDouble())


		###TODO replace this "model" for real one once available in kratos core
		self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

		#construct the solver (main setting methods are located in the solver_module)
		solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
		self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])


		#### Output settings start ####

		self.problem_path = os.getcwd()
		self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()


	def InitializeSolutionStep(self):

		neighbour_elemental_finder =  KratosMultiphysics.FindElementalNeighboursProcess(self.main_model_part,2, 5)
		neighbour_elemental_finder.Execute()

		self.delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

		self.time = self.time + self.delta_time
		self.step = self.step + 1

		self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.step
		self.main_model_part.CloneTimeStep(self.time) 


		print(" [STEP:",self.step," TIME:",self.time,"]")

		# processes to be executed at the begining of the solution step
		self.model_processes.ExecuteInitializeSolutionStep()

		self.GraphicalOutputExecuteInitializeSolutionStep()

		self.solver.InitializeSolutionStep()

	def SolveSolutionStep(self):

		self.clock_time = self.StartTimeMeasuring();

		self.solver.Solve()

		# ************* PRINTS DE PRUEBA ******************************************* 
		print("*************************************")
		print(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
		new_list = []
		for node in self.main_model_part.Nodes:
			pass
			#if node.Id == 3:
				#print("Id: ",node.Id)
			#	print("X: ",node.X)
				#print("Y: ",node.Y)
			#	print("Displ: ", node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT))


		damaged_elements = []
		for elem in self.main_model_part.Elements:
			damage = elem.GetValuesOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.main_model_part.ProcessInfo)
			if damage[0][0] > 0.0:
				print("Id elemento dañado: ", elem.Id)
				Wait()
				#print(damaged_elements)

			#print("damage : ", damage)
			#print("\n")

		print("*************************************")


		#Wait()

		self.StopTimeMeasuring(self.clock_time,"Solving", False);


	def FinalizeSolutionStep(self):


		self.GraphicalOutputExecuteFinalizeSolutionStep()            

		# processes to be executed at the end of the solution step
		self.model_processes.ExecuteFinalizeSolutionStep()

		# processes to be executed before witting the output      
		self.model_processes.ExecuteBeforeOutputStep()

		# write output results GiD: (frequency writing is controlled internally)
		self.GraphicalOutputPrintOutput()            

		# processes to be executed after witting the output
		self.model_processes.ExecuteAfterOutputStep()

		self.main_model_part.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)





# Main de DEM application

class DEM_Solution(MainDEM.Solution):

	def Info(self):
		print("DEM part of the FEM-DEM application")



