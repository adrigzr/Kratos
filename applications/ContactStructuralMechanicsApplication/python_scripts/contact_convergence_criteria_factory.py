## This script collects the available convergence criteria to be used in the SolidMechanicsApplication

from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

try: 
  import KratosMultiphysics.MeshingApplication as MeshingApplication 
  missing_meshing_dependencies = False 
  missing_application = '' 
except ImportError as e: 
    missing_meshing_dependencies = True 
    # extract name of the missing application from the error message 
    import re 
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1) 

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Convergence criteria class
class convergence_criterion:
    def __init__(self, convergence_criterion_parameters, model_part = None, adaptative_parameters = None, processes_list = None):
        # Note that all the convergence settings are introduced via a Kratos parameters object.
        
        if "Contact" in convergence_criterion_parameters["convergence_criterion"].GetString():
            D_RT = convergence_criterion_parameters["displacement_relative_tolerance"].GetDouble()
            D_AT = convergence_criterion_parameters["displacement_absolute_tolerance"].GetDouble()
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            CD_RT = convergence_criterion_parameters["contact_displacement_relative_tolerance"].GetDouble()
            CD_AT = convergence_criterion_parameters["contact_displacement_absolute_tolerance"].GetDouble()
            CR_RT = convergence_criterion_parameters["contact_residual_relative_tolerance"].GetDouble()
            CR_AT = convergence_criterion_parameters["contact_residual_absolute_tolerance"].GetDouble()
            contact_tolerance = convergence_criterion_parameters["contact_tolerance"].GetDouble()
            fancy_convergence_criterion = convergence_criterion_parameters["fancy_convergence_criterion"].GetBool()
            print_convergence_criterion = convergence_criterion_parameters["print_convergence_criterion"].GetBool()
            ensure_contact = convergence_criterion_parameters["ensure_contact"].GetBool()
            echo_level = convergence_criterion_parameters["echo_level"].GetInt()
            
            if(echo_level >= 1):
                print("::[Mechanical Solver]:: CONVERGENCE CRITERION : ", convergence_criterion_parameters["convergence_criterion"].GetString())
            
            if (fancy_convergence_criterion == True):
                table = KratosMultiphysics.TableStreamUtility()
            
            if(convergence_criterion_parameters["convergence_criterion"].GetString() == "Contact_Displacement_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Contact_Residual_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Contact_Mixed_criterion"):
                if (fancy_convergence_criterion == True):
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierMixedContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
                    
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Contact_And_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(convergence_criterion_parameters["convergence_criterion"].GetString() == "Contact_Or_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT, ensure_contact, table, print_convergence_criterion)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact, table, print_convergence_criterion)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, CD_RT, CD_AT,ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, CR_RT, CR_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            
            # Adding the mortar criteria
            if  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictionless"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria(contact_tolerance, table, print_convergence_criterion)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria(contact_tolerance)
            elif  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictional"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria(contact_tolerance, table, print_convergence_criterion)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria(contact_tolerance)
            elif ("MeshTying" in convergence_criterion_parameters["mortar_type"].GetString()):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria(table)
                else:
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
            
            Mortar.SetEchoLevel(echo_level)

            if (fancy_convergence_criterion == True):
                self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar, table, print_convergence_criterion)
            else:
                self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(self.mechanical_convergence_criterion, Mortar)
            self.mechanical_convergence_criterion.SetEchoLevel(echo_level)
            self.mechanical_convergence_criterion.SetActualizeRHSFlag(True)
        
        else: # Standard criteria (same as structural mechanics application)
            # Construction of the class convergence_criterion
            import convergence_criteria_factory
            base_mechanical_convergence_criterion = convergence_criteria_factory.convergence_criterion(convergence_criterion_parameters)
        
            # Adding the mortar criteria
            if  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictionless"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif  (convergence_criterion_parameters["mortar_type"].GetString() == "ALMContactFrictional"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif ("MeshTying" in convergence_criterion_parameters["mortar_type"].GetString()):
                Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
                Mortar.SetEchoLevel(echo_level)
                self.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( base_mechanical_convergence_criterion.mechanical_convergence_criterion, Mortar)
                (self.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
        
        # add remeshing criterion
        #if (convergence_criterion_parameters.Has("remeshing") == True):
        #    Remeshing = MeshingApplication.ErrorMeshCriteria(model_part, adaptative_parameters, processes_list)
        #    self.mechanical_convergence_criterion = ContactStructuralMechanicsApplication.RemeshingAndConvergenceCriteria(self.mechanical_convergence_criterion, Remeshing)

