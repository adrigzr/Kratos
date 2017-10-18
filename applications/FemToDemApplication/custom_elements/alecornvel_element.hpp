//
//   Project Name:        KratosFemToDemApplication $
//   Created by:          $Author:Alejandro Cornejo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                Sept 2016 $
//   Revision:            $Revision:                  0.0 $
//

#if !defined(KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED )
#define  KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED

#include "custom_constitutive/zarate_law.hpp"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "custom_elements/small_displacement_element.hpp"

namespace Kratos
{


	class AleCornVelElement : public SmallDisplacementElement    // Derived Element from SolidMechanics
	{

	public:

		/// Default constructors
		AleCornVelElement(IndexType NewId, GeometryType::Pointer pGeometry);

		AleCornVelElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

		///Copy constructor
		AleCornVelElement(AleCornVelElement const& rOther);

		/// Destructor.
		virtual ~AleCornVelElement();

		/// Assignment operator.
		AleCornVelElement& operator=(AleCornVelElement const& rOther);


		Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


		Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

		AleCornVelElement()
		{
		}
		
		// *************** Methods Alejandro Cornejo ***************
		//**********************************************************

		void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
		void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
		void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);
		void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix, const double &rYoungModulus,
			const double &rPoissonCoefficient);

		void CalculateDN_DX(Matrix& rDN_DX, int PointNumber);

		void CalculateInfinitesimalStrain(Vector& rStrainVector, const Matrix& rDN_DX);

		void CalculateStressVector(Vector& rStressVector, const Matrix& rConstitutiveMAtrix, const Vector& rInfinitesimalStrainVector);

		void CalculatePrincipalStress(Vector& PrincipalStressVector, const Vector StressVector);

		void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

		void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
		void CalculateLocalSystem
			(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo);

		void AverageVector(Vector& rAverageVector, const Vector& v, const Vector& w);

		void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
			std::vector<Vector>& rValues,
			const ProcessInfo& rCurrentProcessInfo);

		void Get2MaxValues(Vector& MaxValues, double a, double b, double c);

		void IntegrateStressDamageMechanics(Vector& rIntegratedStress, 
			double& damage, const Vector StrainVector, const Vector StressVector, int cont, double l_char);


		void ModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char);
		void RankineCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char);
		void DruckerPragerCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char);
		void SimoJuCriterion(Vector& rIntegratedStress, double& damage, const Vector& StrainVector, const Vector& StressVector, int cont, double l_char);

		void TangentModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char);

		// Stress Invariants in 2D
		double Calculate_I1_Invariant(double sigma1, double sigma2);
		double Calculate_J2_Invariant(double sigma1, double sigma2);
		double Calculate_J3_Invariant(double sigma1, double sigma2, double I1);

		void CalculateIntegratedStressVector(Vector& rIntegratedStressVector,const Vector& StressVector, const double damage)
		{
			Vector res = ZeroVector(3);
			res = (1 - damage)*StressVector;
		}

		// Lode's angle
		double Calculate_Theta_Angle(double J2, double J3);

		// Converged values
		void   Set_threshold(double af, int cont) { thresholds[cont] = af; }
		double Get_threshold(int cont) { return thresholds[cont]; }

		void   Set_threshold(double af) { threshold = af; }
		double Get_threshold() { return threshold; }

		void   Set_Convergeddamage(double af) { damage = af; }
		double Get_Convergeddamage() { return damage; }

		void   SetConverged_f_sigma(double af) { f_sigma = af; }
		double GetConverged_f_sigma() { return f_sigma; }

		void   SetConverged_f_sigmas(double af, int cont) { f_sigmas[cont] = af; }
		double GetConverged_f_sigmas(int cont) { return f_sigmas[cont]; }

		void   Set_Convergeddamages(double af, int cont) { damages[cont] = af; }
		double Get_Convergeddamages(int cont) { return damages[cont]; }

		// Non Converged values
		void   Set_NonConvergeddamages(double af, int cont) { NonConvergeddamages[cont] = af; }
		double Get_NonConvergeddamages(int cont) { return NonConvergeddamages[cont]; }

		void   Set_NonConvergeddamage(double af) { NonConvergeddamage = af; }
		double Get_NonConvergeddamage() { return NonConvergeddamage; }

		void   Set_NonConvergedf_sigma(double af, int cont) { NonConverged_f_sigma[cont] = af; }
		double Get_NonConvergedf_sigma(int cont) { return NonConverged_f_sigma[cont]; }

		void   Set_NonConvergedf_sigma(double af) { NonConvergedf_sigma = af; }
		double Get_NonConvergedf_sigma() { return NonConvergedf_sigma; }

		void   ResetNonConvergedVars()
		{
			this->Set_NonConvergeddamage(0.0);
			this->Set_NonConvergedf_sigma(0.0);

			for (int cont = 0;cont < 3;cont++)
			{
				this->Set_NonConvergeddamages(0, cont);
				this->Set_NonConvergedf_sigma(0, cont);
			}
		}

		// Characteristic length Calculations
		void   Set_l_char(double af, int cont) { l_char[cont] = af; }
		double Get_l_char(int cont) { return l_char[cont]; }

		void   SetJ(double af) { Jac = af; }
		double GetJ() { return Jac; }

		// Auxiliar functions...
		void IterationPlus() { iteration++; }
		int  GetIteration() { return iteration; }
		void SetToZeroIteration() { iteration = 0; }
		//void AssignSmoothedStress(Element& Elem);

		void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
		Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

		// Functions to calculate the Constitutive tangent tensor by numerical derivation
		double GetMaxValue(Vector Strain);
		double GetMaxAbsValue(Vector Strain);
		double GetMinAbsValue(Vector Strain);
		void PerturbateStrainComponent(const Vector& rStrainVector, Vector& PertubatedStrain, const double& perturbation, int component);
		double CalculatePerturbation(const Vector& StrainVector, int component);
		void CalculateTangentTensor(Matrix& rTangentTensor, const Vector& StrainVector, const Vector& IntegratedStressVector ,int cont, double l_char);

		void SetStressVector(Vector toStressVector) { toStressVector.resize(3);StressVector = toStressVector; }
		Vector GetStressVector() { return StressVector; }

		void SetStrainVector(Vector toStrainVector) { toStrainVector.resize(3); StrainVector = toStrainVector; }
		Vector GetStrainVector() { return StrainVector; }

		void SetIntegratedStressVector(Vector toIntegratedStressVector) { toIntegratedStressVector.resize(3);IntegratedStressVector = toIntegratedStressVector; }
		Vector GetIntegratedStressVector() { return IntegratedStressVector; }

		void SetBMatrix(Matrix toBMatrix) { toBMatrix.resize(3, 6); B = toBMatrix; }
		Matrix GetBMatrix(){ return B; }

		void CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX);

		void SetIntegrationCoefficient(double tomIntegrationCoefficient){ mIntegrationCoefficient = tomIntegrationCoefficient; }
		double GetIntegrationCoefficient(){ return mIntegrationCoefficient; }
		
	private:
		int iteration = 0;

		// Each component == Each edge
		double f_sigmas[4] = { 0.0, 0.0, 0.0,0.0 };   // Mohr-Coulomb equivalent stress
		double thresholds[4] = { 0.0, 0.0, 0.0, 0.0 };   // Stress Threshold on edge

		double threshold = 0.0;
		double f_sigma   = 0.0;

		double damages[4] = { 0.0, 0.0, 0.0, 0.0 };     // Converged damage on each edge
		double damage = 0.0;                            // Converged damage

		double NonConvergeddamages[4] = { 0.0, 0.0, 0.0, 0.0 };    // Damages on edges of "i" iteration
		double NonConverged_f_sigma[4] = { 0.0, 0.0, 0.0, 0.0 };   // Equivalent stress of "i" iteration

		double NonConvergedf_sigma = 0.0;
		double NonConvergeddamage = 0.0;       // Damage of the element of "i" iteration

		double l_char[3] = { 0.0, 0.0, 0.0 };  // Characteristic length on each edge

		Vector StressVector = ZeroVector(3);
		Vector StrainVector = ZeroVector(3);
		Vector IntegratedStressVector = ZeroVector(3);
		Matrix B = ZeroMatrix(3, 6);

		double Jac = 0.0;
		double mIntegrationCoefficient
	};// Class AleCornVelElement
	
}// Namespace Kratos
#endif // KRATOS_ALECORNVEL_ELEMENT_H_INCLUDED  defined 