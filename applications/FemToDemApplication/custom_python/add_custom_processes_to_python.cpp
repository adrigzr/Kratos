//   
//   Project Name:           
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: $
//#include <boost/python.hpp>

// Project includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "includes/define.h"
#include "processes/process.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_custom_processes_to_python.h"

#include "processes/find_elements_neighbours_process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
//#include "includes/node.h"
//#include "includes/define.h"
//#include "processes/process.h"
// #include "spaces/ublas_space.h"
// #include "linear_solvers/linear_solver.h"
// #include "custom_python/add_custom_processes_to_python.h"


//#include "custom_processes/adaptive_mesh_refinement_process.hpp"
//#include "custom_processes/mapping_variables_process.hpp"  //TODO when p.type done uncomment

namespace Kratos
{

	namespace Python
	{
		
	using namespace boost::python;

	void AddCustomProcessesToPython()
		{
			typedef Process                           ProcessBaseType;
			//typedef AdaptiveMeshRefinementProcess     AdaptiveMeshRefinementProcessType;
			//typedef MappingVariablesProcess           MappingVariablesProcessType;

			class_<FindElementalNeighboursProcess, bases<ProcessBaseType>, boost::noncopyable >
				("FindElementalNeighboursProcess", init<ModelPart&, int, unsigned int>())
				.def("Execute", &FindElementalNeighboursProcess::Execute);
				
			// Adaptive Mesh Refinement Process
			// class_< AdaptiveMeshRefinementProcessType, bases< ProcessBaseType >, boost::noncopyable >
			//    ( "AdaptiveMeshRefinementProcess",
			//    init < ModelPart&,std::string,std::string,std::string,std::string,double,int >());
			  //.def("Execute", &AdaptiveMeshRefinementProcess::Execute);
				
			// // Mapping Variables Process
			// class_< MappingVariablesProcessType, bases< ProcessBaseType >, boost::noncopyable > ( "MappingVariablesProcess",
			// 	init < ModelPart&,ModelPart&,std::string >());


		}

	}  // namespace Python.

} // Namespace Kratos