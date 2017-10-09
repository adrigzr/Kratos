from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
from KratosMultiphysics.EmpireApplication import *
import ctypes as ctp
import os

CheckForPreviousImport()

class EmpireWrapper:
    # Source of Implementation: https://code.activestate.com/recipes/52558/
    # storage for the instance reference
    __instance = None

    def __init__(self):
        """ Create singleton instance """
        # Check whether we already have an instance
        if EmpireWrapper.__instance is None:
            # Create and remember instance
            EmpireWrapper.__instance = EmpireWrapper.__EmpireWrapper()

    def __getattr__(self, attr):
        """ Delegate access to implementation """
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):
        """ Delegate access to implementation """
        return setattr(self.__instance, attr, value)

    class __EmpireWrapper:
        # Wrapper for the EMPIRE API (/EMPIRE-Core/EMPIRE_API/src/include/EMPIRE_API.h)
        # Implemented as Singleton, bcs otherwise the EMPIRE library can be imported several times
        ##### Constructor #####
        # -------------------------------------------------------------------------------------------------
        def __init__(self):
            self.model_parts = {}
            self._LoadEmpireLibrary()
        # -------------------------------------------------------------------------------------------------

        ##### Public Functions #####
        # -------------------------------------------------------------------------------------------------
        def Connect(self, xml_input_file):
            ''' Establishes the necessary connection with the Emperor '''
            self.libempire_api.EMPIRE_API_Connect(xml_input_file.encode())
            print("::EMPIRE:: Connected")
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def Disconnect(self):
            ''' Performs disconnection and finalization operations to the Emperor '''
            self.libempire_api.EMPIRE_API_Disconnect()
            print("::EMPIRE:: Disconnected")
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendMesh(self, mesh_name, model_part):
            ''' Send the mesh to the server
            \param[in] name name of the mesh
            \param[in] numNodes number of nodes
            \param[in] numElems number of elements
            \param[in] nodes coordinates of all nodes
            \param[in] nodeIDs IDs of all nodes
            \param[in] numNodesPerElem number of nodes per element
            \param[in] elems connectivity table of all elements
            
            void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs,
                    int *numNodesPerElem, int *elems); '''
            # mesh_name: name of mesh in the emperor input
            
            # Save the ModelPart for data-field exchange later
            self._SaveModelPart(mesh_name, model_part)
            
            # extract interface mesh information
            numNodes = [];          numElems = []
            nodeCoors = [];         nodeIDs = []
            numNodesPerElem = [];   elemTable = []
            self._GetMesh(model_part, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)

            # convert python lists to ctypes, required for empire-function call
            c_numNodes = (ctp.c_int * len(numNodes))(*numNodes)
            c_numElems = (ctp.c_int * len(numElems))(*numElems)
            c_nodeCoors = (ctp.c_double * len(nodeCoors))(*nodeCoors)
            c_nodeIDs = (ctp.c_int * len(nodeIDs))(*nodeIDs)
            c_numNodesPerElem = (ctp.c_int * len(numNodesPerElem))(*numNodesPerElem)
            c_elemTable = (ctp.c_int * len(elemTable))(*elemTable)

            # send mesh to Emperor
            self.libempire_api.EMPIRE_API_sendMesh("mesh_name.encode()", 
                                                c_numNodes[0], c_numElems[0], 
                                                c_nodeCoors, c_nodeIDs, 
                                                c_numNodesPerElem, c_elemTable)
            print("::EMPIRE:: Sent Mesh")
        # -------------------------------------------------------------------------------------------------
        
        # -------------------------------------------------------------------------------------------------
        def ReceiveMesh(self, mesh_name, model_part):
            ''' Recieve mesh from the server
            \param[in] name name of the mesh
            \param[in] numNodes number of nodes
            \param[in] numElems number of elements
            \param[in] nodes coordinates of all nodes
            \param[in] nodeIDs IDs of all nodes
            \param[in] numNodesPerElem number of nodes per element
            \param[in] elems connectivity table of all elements
            
            void EMPIRE_API_recvMesh(char *name, int *numNodes, int *numElems, double **nodes, int **nodeIDs,
                    int **numNodesPerElem, int **elem); '''
            # mesh_name: name of mesh in the emperor input

            # Save the ModelPart for data-field exchange later
            self._SaveModelPart(mesh_name, model_part)

            c_numNodes = ctp.pointer(ctp.c_int(0))
            c_numElems = ctp.pointer(ctp.c_int(0))
            c_nodeCoors = ctp.pointer(ctp.pointer(ctp.c_double(0)))
            c_nodeIDs = ctp.pointer(ctp.pointer(ctp.c_int(0)))
            c_numNodesPerElem = ctp.pointer(ctp.pointer(ctp.c_int(0)))
            c_elemTable = ctp.pointer(ctp.pointer(ctp.c_int(0)))
            
            self.libempire_api.EMPIRE_API_recvMesh(mesh_name.encode(), 
                                                   c_numNodes, c_numElems, 
                                                   c_nodeCoors, c_nodeIDs, 
                                                   c_numNodesPerElem, c_elemTable)

            numNodes = c_numNodes.contents.value
            numElems = c_numElems.contents.value
            nodeCoors = c_nodeCoors.contents
            nodeIDs = c_nodeIDs.contents
            numNodesPerElem = c_numNodesPerElem.contents
            elemTable = c_elemTable.contents
            
            self._SetMesh(model_part, numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable)
            print("::EMPIRE:: Received Mesh")
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def SendDataField(self, mesh_name, data_field_name, kratos_variable):
            ''' Send data field to the server
            \param[in] name name of the field
            \param[in] sizeOfArray size of the array (data field)
            \param[in] dataField the data field to be sent
            
            void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField); '''
            # mesh_name: name of mesh in the emperor input
            # data_field_name: name of dataField in the emperor input

            # get ModelPart
            model_part = self.model_parts[mesh_name]

            # extract data field from nodes
            data_field = []
            self._GetDataField(model_part, kratos_variable, data_field)

            # convert list containg the data field to ctypes
            c_data_field = (ctp.c_double * len(data_field))(*data_field)
            c_size = len(c_data_field)

            # send data field to EMPIRE
            self.libempire_api.EMPIRE_API_sendDataField("data_field_name.encode()", c_size, c_data_field)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveDataField(self, mesh_name, data_field_name, kratos_variable):   
            ''' Receive data field from the server
            \param[in] name name of the field
            \param[in] sizeOfArray size of the array (data field)
            \param[out] dataField the data field to be received
            
            void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField); '''
            # mesh_name: name of mesh in the emperor input
            # data_field_name: name of dataField in the emperor input

            # get ModelPart
            model_part = self.model_parts[mesh_name]

            # Determine Size of Variable
            size_of_variable = self._SizeOfVariable(model_part, kratos_variable)

            # initialize vector storing the values
            size_data_field = model_part.NumberOfNodes() * size_of_variable
            c_size_data_field = ctp.c_int(size_data_field)
            c_data_field = (ctp.c_double * size_data_field)(0)

            # receive data field from empire
            self.libempire_api.EMPIRE_API_recvDataField(data_field_name.encode(), c_size_data_field, c_data_field)

            self._SetDataField(model_part, kratos_variable, c_data_field, size_of_variable)
        # -------------------------------------------------------------------------------------------------
        
        # -------------------------------------------------------------------------------------------------
        def SendArray(self, array_name, array_to_send):
            ''' Send signal to the server
            \param[in] name name of the signal
            \param[in] sizeOfArray size of the array (signal)
            \param[in] signal the signal
            
            void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal); '''
            # array_name: name of signal in the emperor input

            # convert array to ctypes
            c_signal = (ctp.c_double * len(array_to_send))(*array_to_send)
            c_size = len(c_signal)

            self.libempire_api.EMPIRE_API_sendSignal_double(array_name.encode(), c_size, c_signal)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveArray(self, array_name, array_size):
            ''' Receive signal from the server
            \param[in] name name of the signal
            \param[in] sizeOfArray size of the array (signal)
            \param[in] signal the signal
            
            void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal); '''
            # array_name: name of signal in the emperor input

            # initialize vector storing the values
            c_signal = (ctp.c_double * array_size)(0)

            self.libempire_api.EMPIRE_API_recvSignal_double(array_name.encode(), array_size, c_signal)

            return self._ConvertToList(array_size, c_signal)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def ReceiveConvergenceSignal(self):
            '''Receive the convergence signal of an loop
            \return 1 means convergence, 0 means non-convergence
            
            int EMPIRE_API_recvConvergenceSignal(); '''

            return self.libempire_api.EMPIRE_API_recvConvergenceSignal()
        # -------------------------------------------------------------------------------------------------

        ##### Private Functions #####
        # -------------------------------------------------------------------------------------------------
        def _LoadEmpireLibrary(self):
            if hasattr(self, 'libempire_api'): # the library has been loaded already
                raise ImportError("The EMPIRE library must be loaded only once!")

            if "EMPIRE_API_LIBSO_ON_MACHINE" not in os.environ:
                raise ImportError("The EMPIRE environment is not set!")

            try: # OpenMPI
                self.libempire_api = ctp.CDLL(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'], ctp.RTLD_GLOBAL)
                print("::EMPIRE:: Using standard OpenMPI")
            except: # Intel MPI or OpenMPI compiled with "–disable-dlopen" option
                self.libempire_api = ctp.cdll.LoadLibrary(os.environ['EMPIRE_API_LIBSO_ON_MACHINE'])
                print("::EMPIRE:: Using Intel MPI or OpenMPI compiled with \"–disable-dlopen\" option")
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _GetMesh(self, model_part, num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table):
            num_nodes.append(model_part.NumberOfNodes())
            num_elements.append(model_part.NumberOfElements())
            
            for node in model_part.Nodes:
                node_coords.append(node.X)
                node_coords.append(node.Y)
                node_coords.append(node.Z)
                node_IDs.append(node.Id)
                
            for elem in model_part.Elements:
                num_nodes_per_element.append(len(elem.GetNodes()))
                for node in elem.GetNodes():
                    element_table.append(node.Id)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _SetMesh(self, model_part, num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table):
            # This function requires an empty ModelPart
            # It constructs Nodes and Elements from what was received from EMPIRE

            # Some checks to validate input:
            if model_part.NumberOfNodes() != 0:
                raise Exception("ModelPart is not empty, it has some Nodes!")
            if model_part.NumberOfElements() != 0:
                raise Exception("ModelPart is not empty, it has some Elements!")
            if model_part.NumberOfConditions() != 0:
                raise Exception("ModelPart is not empty, it has some Conditions!")

            # Create Nodes
            for i in range(num_nodes):
                model_part.CreateNewNode(node_IDs[i], node_coords[3*i+0], node_coords[3*i+1], node_coords[3*i+2]) # Id, X, Y, Z

            # Create dummy Property for Element
            model_part.AddProperties(Properties(1))
            prop = model_part.GetProperties()[1]

            element_table_counter = 0
            # Create Elements
            for i in range(num_elements):
                num_nodes_element = num_nodes_per_element[i]
                if num_nodes_element == 2:
                    name_element = "Element2D2N"
                elif num_nodes_element == 3:
                    name_element = "Element2D3N"
                elif num_nodes_element == 4: # TODO how to distinguish from Tetras?
                    name_element = "Element2D4N"
                else:
                    raise Exception("Wrong number of nodes for creating the element")

                element_nodes = []
                for j in range(num_nodes_element):
                    element_nodes.append(int(element_table[element_table_counter]))
                    element_table_counter += 1
                    
                model_part.CreateNewElement(name_element, i+1, element_nodes, prop)
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _GetDataField(self, model_part, kratos_variable, data_field):
            size_of_variable = self._SizeOfVariable(model_part, kratos_variable)

            for node in model_part.Nodes:
                data_value = node.GetSolutionStepValue(kratos_variable)

                if size_of_variable == 1:
                    data_field.append(data_value)
                else:
                    for i in range(size_of_variable):
                        data_field.append(data_value[i])
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _SetDataField(self, model_part, kratos_variable, data_field, size_of_variable):
            # check if size of data field is correct
            if len(data_field) != model_part.NumberOfNodes() * size_of_variable:
                raise("ERROR: received data field has wrong size!")

            if size_of_variable > 1:
                value = Vector(size_of_variable)

            i = 0
            # assign values to nodes of interface for current time step
            for node in model_part.Nodes:
                if size_of_variable == 1:
                    value = data_field[size_of_variable * i]
                else:
                    for j in range(size_of_variable):
                        value[j] = data_field[size_of_variable * i + j]
                
                node.SetSolutionStepValue(kratos_variable, 0, value)
                
                i = i + 1
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _SizeOfVariable(self, model_part, kratos_variable):
            # this function is very general, even though EMPIRE works with Scalar and Vector quantities only!
            try:
                first_node = next(iter(model_part.Nodes))
                value = first_node.GetSolutionStepValue(kratos_variable)
                if (isinstance(value, float) or isinstance(value, int)): # Variable is a scalar
                    size_of_variable = 1
                else:
                    size_of_variable = len(first_node.GetSolutionStepValue(kratos_variable))
            except StopIteration:
                raise TypeError("size_of_variable could not be determined")
            
            return size_of_variable
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _SaveModelPart(self, mesh_name, model_part):
            # Save the model_part for data-field exchange later
            if mesh_name in self.model_parts:
                raise ValueError("Mesh exsts already")
            else:
                self.model_parts.update({mesh_name : model_part})
        # -------------------------------------------------------------------------------------------------

        # -------------------------------------------------------------------------------------------------
        def _ConvertToList(self, array_size, c_signal):
            converted_list = [0.0] * array_size # preallocate

            for i in range(array_size):
                converted_list[i] = c_signal[i]
            
            return converted_list
        # -------------------------------------------------------------------------------------------------
