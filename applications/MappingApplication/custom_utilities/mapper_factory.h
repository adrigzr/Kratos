//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_FACTORY_H_INCLUDED )
#define  KRATOS_MAPPER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/nearest_neighbor_mapper.h"
#include "custom_utilities/nearest_element_mapper.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Python Interface of the MappingApplication
/** This class constructs the mappers and exposes them to Python
* Some checks are performed to see if the Input (ModelParts and JSON-Parameters) are valid
* Also the additional timing information is implemented here (echo_level = 1)
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
class MapperFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperFactory
    KRATOS_CLASS_POINTER_DEFINITION(MapperFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperFactory(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                  Parameters JsonParameters) :
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mJsonParameters(JsonParameters)
    {
        ReadInterfaceModelParts();
        ConstructMapper();
    }

    /// Destructor.
    virtual ~MapperFactory() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags& rOptions, double SearchRadius)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->UpdateInterface(rOptions, SearchRadius);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "UpdateInterface",
                elapsed_time);
    }


    /* This function maps a variable from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "Map",
                elapsed_time);
    }

    /* This function maps a variable from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->Map(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "Map",
                elapsed_time);
    }


    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "InverseMap",
                elapsed_time);
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags& rOptions)
    {
        double start_time = MapperUtilities::GetCurrentTime();
        mpMapper->InverseMap(rOriginVariable, rDestinationVariable, rOptions);
        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "InverseMap",
                elapsed_time);
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MapperFactory" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperFactory";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Mapper::Pointer mpMapper;
    std::string mMapperType;

    ModelPart& mrModelPartOrigin;
    ModelPart& mrModelPartDestination;

    ModelPart* mpInterfaceModelPartOrigin;
    ModelPart* mpInterfaceModelPartDestination;

    Parameters mJsonParameters;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ReadInterfaceModelParts()
    {
        int echo_level = 0;
        // read the echo_level temporarily, bcs the mJsonParameters have not yet been validated and defaults assigned
        if (mJsonParameters.Has("echo_level"))
        {
            echo_level = std::max(echo_level, mJsonParameters["echo_level"].GetInt());
        }

        int comm_rank_origin = mrModelPartOrigin.GetCommunicator().MyPID();
        int comm_rank_destination = mrModelPartDestination.GetCommunicator().MyPID();

        if (mJsonParameters.Has("interface_submodel_part_origin"))
        {
            std::string name_interface_submodel_part = mJsonParameters["interface_submodel_part_origin"].GetString();
            mpInterfaceModelPartOrigin = &mrModelPartOrigin.GetSubModelPart(name_interface_submodel_part);

            if (echo_level >= 3 && comm_rank_origin == 0)
            {
                std::cout << "Mapper: SubModelPart used for Origin-ModelPart" << std::endl;
            }
        }
        else
        {
            mpInterfaceModelPartOrigin = &mrModelPartOrigin;

            if (echo_level >= 3 && comm_rank_origin == 0)
            {
                std::cout << "Mapper: Main ModelPart used for Origin-ModelPart" << std::endl;
            }
        }

        if (mJsonParameters.Has("interface_submodel_part_destination"))
        {
            std::string name_interface_submodel_part = mJsonParameters["interface_submodel_part_destination"].GetString();
            mpInterfaceModelPartDestination = &mrModelPartDestination.GetSubModelPart(name_interface_submodel_part);

            if (echo_level >= 3 && comm_rank_destination == 0)
            {
                std::cout << "Mapper: SubModelPart used for Destination-ModelPart" << std::endl;
            }
        }
        else
        {
            mpInterfaceModelPartDestination = &mrModelPartDestination;

            if (echo_level >= 3 && comm_rank_destination == 0)
            {
                std::cout << "Mapper: Main ModelPart used for Destination-ModelPart" << std::endl;
            }
        }
    }

    void ConstructMapper()
    {
        double start_time = MapperUtilities::GetCurrentTime();

        if (!mJsonParameters.Has("mapper_type"))
        {
            KRATOS_ERROR << "No \"mapper_type\" defined in json" << std::endl;
        }

        mMapperType = mJsonParameters["mapper_type"].GetString();

        if (mMapperType == "NearestNeighbor")
        {
            if (mJsonParameters.Has("approximation_tolerance"))
            {
                KRATOS_ERROR << "Invalid Parameter \"approximation_tolerance\" "
                             << "specified for Nearest Neighbor Mapper" << std::endl;
            }

            mpMapper = Mapper::Pointer(new NearestNeighborMapper(*mpInterfaceModelPartOrigin,
                                       *mpInterfaceModelPartDestination,
                                       mJsonParameters));
        }
        else if (mMapperType == "NearestElement")
        {
            mpMapper = Mapper::Pointer(new NearestElementMapper(*mpInterfaceModelPartOrigin,
                                       *mpInterfaceModelPartDestination,
                                       mJsonParameters));

        } /*else if (mMapperType == "Barycentric") {
              mpMapper = Mapper::Pointer(new BarycentricMapper(*mpInterfaceModelPartOrigin,
                                                                 *mpInterfaceModelPartDestination,
                                                                 mJsonParameters));

          } *//*else if (mMapperType == "RBF") {
              mpMapper = Mapper::Pointer(new RBFMapper(*mpInterfaceModelPartOrigin,
                                                         *mpInterfaceModelPartDestination,
                                                         mJsonParameters));

          } *//*else if (mMapperType == "Mortar") {
              mpMapper = Mapper::Pointer(new MortarMapper(*mpInterfaceModelPartOrigin,
                                                            *mpInterfaceModelPartDestination,
                                                            mJsonParameters));

          } *//*else if (mMapperType == "IGA") {
              mpMapper = Mapper::Pointer(new IGAMapper(*mpInterfaceModelPartOrigin,
                                                         *mpInterfaceModelPartDestination,
                                                         mJsonParameters));

          } */else
        {
            KRATOS_ERROR << "Selected Mapper \"" << mMapperType << "\" not implemented" << std::endl;
        }

        double elapsed_time = MapperUtilities::GetCurrentTime() - start_time;

        mpMapper->pGetMapperCommunicator()->PrintTime(mMapperType,
                "Mapper Construction",
                elapsed_time);
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{s


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MapperFactory& operator=(MapperFactory const& rOther);

    //   /// Copy constructor.
    //   MapperFactory(MapperFactory const& rOther){}


    ///@}

}; // Class MapperFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapperFactory& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapperFactory& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_FACTORY_H_INCLUDED  defined