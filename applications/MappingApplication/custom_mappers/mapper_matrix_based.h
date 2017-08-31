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

#if !defined(KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED )
#define  KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_utilities/interface_preprocess.h"
// #include "custom_strategies/builders/mapping_matrix_builder.h"
// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
// #include "custom_strategies/builders/trilinos_mapping_matrix_builder.h"
// #endif

namespace Kratos
{

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
/// Short class definition.
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TMappingMatrixBuilder,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>
class MapperMatrixBased : public Mapper
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MapperMatrixBased
    KRATOS_CLASS_POINTER_DEFINITION(MapperMatrixBased);

    typedef TMappingMatrixBuilder TMappingMatrixBuilderType;

    typedef typename TMappingMatrixBuilderType::Pointer TMappingMatrixBuilderPointerType;

    // typedef typename TMappingMatrixBuilderType::TSparseSpace TSomeOtherSpaceSpace;

    typedef std::unordered_map<int, Node<3> *> EquationIdMapType;

    typedef typename TSparseSpace::DataType TDataType;

    typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;

    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    MapperMatrixBased(ModelPart &rModelPartOrigin, ModelPart &rModelpartDestination,
                      Parameters rJsonParameters) : Mapper(rModelPartOrigin, rModelpartDestination, rJsonParameters)
    {
        mpMdo = TSparseSpace::CreateEmptyMatrixPointer();
        mpQo = TSparseSpace::CreateEmptyVectorPointer();
        mpQd = TSparseSpace::CreateEmptyVectorPointer();


        mpMappingMatrixBuilder = TMappingMatrixBuilderPointerType(new TMappingMatrixBuilderType());


        mpInterfacePreprocessor = InterfacePreprocess::Pointer( new InterfacePreprocess(this->mrModelPartDestination, 
            this->mpInterfaceModelPart) );
    }

    /// Destructor.
    virtual ~MapperMatrixBased() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void UpdateInterfaceSpecific(Kratos::Flags MappingOptions) override {
        if (MappingOptions.Is(MapperFlags::REMESHED))
        {
            ComputeInterfaceModelPart();

            // TODO this gives compiler errors!
            // TSparseSpace::Clear(mpMdo);
            // TSparseSpace::Clear(mpQo);
            // TSparseSpace::Clear(mpQd);
        }
        else
        {
            // TODO this gives compiler errors!
            // TSparseSpace::ClearData(mpMdo);
            // TSparseSpace::ClearData(mpQo);
            // TSparseSpace::ClearData(mpQd);
        }
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {

    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        
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
    virtual std::string Info() const override
    {
        return "MapperMatrixBased";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperMatrixBased";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

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

    TMappingMatrixBuilderPointerType mpMappingMatrixBuilder;

    TSystemVectorPointerType mpQo;
    TSystemVectorPointerType mpQd;
    TSystemMatrixPointerType mpMdo;

    ModelPart::Pointer mpInterfaceModelPart;
    InterfacePreprocess::Pointer mpInterfacePreprocessor;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ComputeMappingMatrix()
    {
        mpMappingMatrixBuilder->SetUpSystem(mrModelPartOrigin, mEquationIdNodeMapOrigin);
        mpMappingMatrixBuilder->SetUpSystem(mrModelPartDestination, mEquationIdNodeMapDestination);

        // mpMappingMatrixBuilder->BuildLHS(BaseType::GetModelPart(), mpMdo);
    }

    /**
    This function creates the Interface-ModelPart and sets up the Matrix-Structure
    */
    virtual void Initialize() = 0;

    /**
    This function creates the Interface-ModelPart
    */
    virtual void ComputeInterfaceModelPart() = 0;

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

    EquationIdMapType mEquationIdNodeMapOrigin; // This is the equivalent to the dofset I think
    EquationIdMapType mEquationIdNodeMapDestination;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{





//     void InitializeMapperStrategy()
//     {
// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
//         if (mCommSize > 1)
//         {
//             mpMapperStrategy = MapperStrategy<TrilinosSparseSpaceType, 
//             LocalSpaceType, TrilinosLinearSolverType>(mrModelPartOrigin, 
//                             mrModelPartDestination);
//         }
//         else
//         {
//             mpMapperStrategy = MapperStrategy<SerialSparseSpaceType, 
//             LocalSpaceType, SerialLinearSolverType>(mrModelPartOrigin, 
//                             mrModelPartDestination);
//         }

// #else // serial compilation
//         mpMapperStrategy = MapperStrategy<SerialSparseSpaceType, 
//         LocalSpaceType, SerialLinearSolverType>(mrModelPartOrigin, 
//                         mrModelPartDestination);
// #endif


// }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    //MapperMatrixBased(MapperMatrixBased const& rOther);

    ///@}

}; // Class MapperMatrixBased

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


// inline std::istream & operator >>(std::istream& rIStream,
//                                   MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     return rIStream;
// }

// /// output stream function

// inline std::ostream & operator <<(std::ostream& rOStream,
//                                   const MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED  defined