//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION)
#define KRATOS_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// An objective function for flow mal distribution.
template <unsigned int TDim>
class MalDistributionObjectiveFunction : public ObjectiveFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MalDistributionObjectiveFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MalDistributionObjectiveFunction(Parameters& rParameters)
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "objective_type": "flow_mal_distribution",
            "surface_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "flow_mal_distribution_direction": "surface_normal"
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mSurfaceModelPartName = rParameters["surface_model_part_name"].GetString();

        if (rParameters["flow_mal_distribution_direction"].IsArray() == false)
            if (rParameters["flow_mal_distribution_direction"].GetString().compare("surface_normal")==0)
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL;
            else
                KRATOS_THROW_ERROR(std::runtime_error,
                               "flow_mal_distribution_direction only accepts as text inputs \"surface_normal\" or \"absolute_value\":",
                               rParameters.PrettyPrintJsonString())
        else
            if (rParameters["flow_mal_distribution_direction"].size() != 3)
            {
                KRATOS_THROW_ERROR(std::runtime_error,
                                   "flow_mal_distribution_direction vector does not have size 3:",
                                   rParameters.PrettyPrintJsonString())
            }
            else
            {
                mMalDistributionDirection = MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION;
                for (unsigned int d = 0; d < TDim; ++d)
                    mDirection[d] = rParameters["flow_mal_distribution_direction"][d].GetDouble();

                if (std::abs(norm_2(mDirection) - 1.0) > 1e-3)
                {
                    const double magnitude = norm_2(mDirection);
                    if (magnitude == 0.0)
                        KRATOS_THROW_ERROR(std::runtime_error,
                                        "flow_mal_distribution_direction is not properly defined.",
                                        "")

                    std::cout
                        << "WARNING: non unit vector detected in \"flow_mal_distribution_direction\": "
                        << rParameters.PrettyPrintJsonString() << std::endl;
                    std::cout << "normalizing \"flow_mal_distribution_direction\"..." << std::endl;

                    for (unsigned int d = 0; d < TDim; d++)
                        mDirection[d] /= magnitude;
                }
            }        

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~MalDistributionObjectiveFunction()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        if (rModelPart.HasSubModelPart(mSurfaceModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"surface_model_part_name\": ",
                mSurfaceModelPartName)
        }


// initialize the variables to zero.
#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
                it->Set(OBJECTIVE_SURFACE, false);
        }

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        // mark objective surface
        for (auto it = rSurfaceModelPart.NodesBegin();
             it != rSurfaceModelPart.NodesEnd();
             ++it)
            it->Set(OBJECTIVE_SURFACE, true);

        KRATOS_CATCH("")
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // allocate auxiliary memory. this is done here instead of Initialize()
        // in case of restart.
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);
        auto normal_calculation_obj = NormalCalculationUtils();
        normal_calculation_obj.CalculateOnSimplex(rSurfaceModelPart, 3);

        CalculateAverageObjectiveParameters(rModelPart);
        mObjectiveValue = Calculate(rModelPart);

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointVelocityContribution(const Element& rElem,
                                                      const Matrix& rAdjointMatrix,
                                                      Vector& rRHSContribution,
                                                      ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        //Question: about the nodes being accounted by this element and the adjacent element as well??

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        const unsigned int NumNodes = rElem.GetGeometry().PointsNumber();

        double inv_coeff = 1/(mAverageSurfaceArea*(mN-1)*mObjectiveValue);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(OBJECTIVE_SURFACE))
                {
                    const array_1d<double,3>& normal = rElem.GetGeometry()[iNode].GetValue(NORMAL);
                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].GetValue(VELOCITY);
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = inv_coeff *
                                                        normal[d] * (
                                                                normal[d]*
                                                                velocity[d] / 
                                                                mAverageSurfaceArea - mAverageVelocity
                                                                );
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            unsigned int LocalIndex = 0;
            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                if (rElem.GetGeometry()[iNode].Is(OBJECTIVE_SURFACE))
                {
                    const array_1d<double,3>& normal = rElem.GetGeometry()[iNode].GetValue(NORMAL);
                    const array_1d<double,3>& velocity = rElem.GetGeometry()[iNode].GetValue(VELOCITY);
                    double magnitude = norm_2(normal);

                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = inv_coeff*
                                                        mDirection[d]*
                                                        magnitude*(
                                                                        mDirection[d] *
                                                                        magnitude *
                                                                        velocity[d] / 
                                                                        mAverageSurfaceArea - mAverageVelocity
                                                                    );
                }
                else
                {
                    for (unsigned int d = 0; d < TDim; d++)
                        rRHSContribution[LocalIndex++] = 0.0;
                }

                //Adding pressure DOF                        
                rRHSContribution[LocalIndex++] = 0.0;
            }
        }

        KRATOS_CATCH("")
    }

    virtual void CalculateAdjointAccelerationContribution(const Element& rElem,
                                                          const Matrix& rAdjointMatrix,
                                                          Vector& rRHSContribution,
                                                          ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rAdjointMatrix.size1())
            rRHSContribution.resize(rAdjointMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual void CalculateSensitivityContribution(const Element& rElem,
                                                  const Matrix& rDerivativesMatrix,
                                                  Vector& rRHSContribution,
                                                  ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rRHSContribution.size() != rDerivativesMatrix.size1())
            rRHSContribution.resize(rDerivativesMatrix.size1(), false);

        for (unsigned int k = 0; k < rRHSContribution.size(); ++k)
            rRHSContribution[k] = 0.0;

        KRATOS_CATCH("")
    }

    virtual double Calculate(ModelPart& rModelPart)
    {
        double result = 0.0;

        CalculateAverageObjectiveParameters(rModelPart);

        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& normal = it->GetValue(NORMAL);
                    const array_1d<double,3>& velocity = it->GetValue(VELOCITY);

                    result += pow(
                                    (
                                        normal[0]*velocity[0] + 
                                        normal[1]*velocity[1] + 
                                        normal[2]*velocity[2]
                                    ) / mAverageSurfaceArea - mAverageVelocity,
                                    2
                                );
                }
            }

        } 
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            #pragma omp parallel reduction(+:result)
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rSurfaceModelPart.Nodes(), NodesBegin, NodesEnd);

                for (auto it = NodesBegin; it != NodesEnd; ++it)
                {
                    const array_1d<double,3>& normal = it->GetValue(NORMAL);
                    const array_1d<double,3>& velocity = it->GetValue(VELOCITY);
                    
                    double magnitude = norm_2(normal);

                    result += pow(
                                    (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    )*magnitude / mAverageSurfaceArea - mAverageVelocity,
                                    2
                                );
                }
            }

        }

        return sqrt(result/(mN-1));
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mSurfaceModelPartName;
    array_1d<double, TDim> mDirection;
    std::vector<Vector> mObjectiveFlagVector;
    //std::vector<unsigned int> mElementIds; //Why do you need it?
    int mMalDistributionDirection;
    double mAverageSurfaceArea;
    double mN;
    double mAverageVelocity;
    double mObjectiveValue;

    const int MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL = 1;
    const int MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION = 2;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    void CalculateAverageObjectiveParameters(ModelPart& rModelPart)
    {
        ModelPart& rSurfaceModelPart = rModelPart.GetSubModelPart(mSurfaceModelPartName);

        mAverageSurfaceArea = 0.0;
        mAverageVelocity = 0.0;
        mN = 0;

        if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_SURFACE_NORMAL)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& normal = it->GetValue(NORMAL);
                double magnitude = norm_2(normal);
                const array_1d<double,3>& velocity = it->GetValue(VELOCITY);

                mAverageVelocity += (
                                        normal[0]*velocity[0] + 
                                        normal[1]*velocity[1] + 
                                        normal[2]*velocity[2]
                                    );
                mAverageSurfaceArea += magnitude;
                mN++;
            }
        }
        else if (mMalDistributionDirection == MAL_DISTRIBUTION_DIRECTION_CUSTOM_DIRECTION)
        {
            for (auto it = rSurfaceModelPart.NodesBegin();
                it != rSurfaceModelPart.NodesEnd();
                ++it)
            {
                const array_1d<double,3>& normal = it->GetValue(NORMAL);
                const array_1d<double,3>& velocity = it->GetValue(VELOCITY);

                double magnitude = norm_2(normal);

                mAverageVelocity += (
                                        mDirection[0]*velocity[0] + 
                                        mDirection[1]*velocity[1] + 
                                        mDirection[2]*velocity[2]
                                    )*magnitude;

                mAverageSurfaceArea += magnitude;
                mN++;
            }
        }

        //Getting the average velocity by dividing flow rate by total objective surface area
        mAverageVelocity /= mAverageSurfaceArea;
        mAverageSurfaceArea /= mN;
    }
    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_FLOW_MAL_DISTRIBUTION_OBJECTIVE_FUNCTION defined */
