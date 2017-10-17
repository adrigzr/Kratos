//
//   Project Name:        KratosDamageApplication $
//   Last modified by:    $Author: Ignasi de Pouplana Sardà $
//   Date:                $Date:                  July 2014 $
//   Revision:            $Revision:                    0.0 $
//

#if !defined(KRATOS_MAPPING_VARIABLES_PROCESS )
#define  KRATOS_MAPPING_VARIABLES_PROCESS

#include <cmath>

#include "includes/model_part.h"
#include "processes/process.h"

#include "fem_to_dem_application_variables.h"

namespace Kratos
{

//Only for Triangles2D3N or Quadrilateral2D4N

class MappingVariablesProcess : public Process
{

protected:

    struct NodeNew
    {
      Element::GeometryType::PointType& rNode;
      int Row,Column;
       
      NodeNew(Element::GeometryType::PointType& NodeReference,int Row_i,int Column_j) : rNode(NodeReference)
      {
        Row = Row_i;
        Column = Column_j;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct ElementOldCell
    {
      std::vector<Element::Pointer> ElementOldVector;
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointNew
    {
      ConstitutiveLaw::Pointer pConstitutiveLaw;
      double X_coord,Y_coord;
      int Row,Column;
      
      GaussPointNew(ConstitutiveLaw::Pointer ConstitutiveLawPointer,double X,double Y,int Row_i,int Column_j)
      {
        pConstitutiveLaw = ConstitutiveLawPointer;
        X_coord = X;
        Y_coord = Y;
        Row = Row_i;
        Column = Column_j;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointOld
    {
      ConstitutiveLaw::Pointer pConstitutiveLaw;
      double X_coord,Y_coord;
      
      GaussPointOld(ConstitutiveLaw::Pointer ConstitutiveLawPointer,double X,double Y)
      {
        pConstitutiveLaw = ConstitutiveLawPointer;
        X_coord = X;
        Y_coord = Y;
      }
    };
    
    //------------------------------------------------------------------------------------
    
    struct GaussPointOldCell
    {
      std::vector<GaussPointOld> GaussPointOldVector;
    };
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

public:
  
    typedef ModelPart::ElementsContainerType ElementsArrayType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Constructor
    MappingVariablesProcess(ModelPart& r_model_part_old, ModelPart& r_model_part_new, std::string imposed_displacement) : mmodel_part_old(r_model_part_old), mmodel_part_new(r_model_part_new)
    {
        mImposedDisplacement = imposed_displacement;
    }
    
    //------------------------------------------------------------------------------------
    
    // Destructor
    virtual ~MappingVariablesProcess(){}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Main Function
    void Execute()
    {
        double X_max = 0.0, X_min = 0.0, Y_max = 0.0, Y_min = 0.0;
    
        double DamageThreshold, CharacteristicLength = 0.0;
    
        this->Initialize(X_max,X_min,Y_max,Y_min,CharacteristicLength,DamageThreshold);
    
        this->NodalDisplacementsMapping(X_max,X_min,Y_max,Y_min);
    
        this->GaussPointStateVariableMapping(X_max,X_min,Y_max,Y_min,CharacteristicLength,DamageThreshold);
    
        this->TransferProcessInfoVariables();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    ModelPart& mmodel_part_old;
    ModelPart& mmodel_part_new;
    
    std::string mImposedDisplacement;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize(double& rX_max,double& rX_min,double& rY_max,double& rY_min, double& rCharacteristicLength, double& rDamageThreshold)
    {
        rX_max = mmodel_part_old.NodesBegin()->X0();
        rX_min = rX_max;
        rY_max = mmodel_part_old.NodesBegin()->Y0();
        rY_min = rY_max;

        double X_me, Y_me;

        for(ModelPart::NodeIterator i = mmodel_part_old.NodesBegin(); i != mmodel_part_old.NodesEnd(); ++i)
        {
            X_me = i->X0();
            Y_me = i->Y0();

            if( X_me > rX_max )  rX_max = X_me;
            else if( X_me < rX_min )  rX_min = X_me;

            if( Y_me > rY_max )  rY_max = Y_me;
            else if( Y_me < rY_min )  rY_min = Y_me;

            //Move old mesh to the original position (to work with both meshes in the reference state)
            (i)->X() = (i)->X0();
            (i)->Y() = (i)->Y0();
            (i)->Z() = (i)->Z0();
        }
        
        rCharacteristicLength = (*(mmodel_part_old.Elements().ptr_begin()))->GetProperties()[CHARACTERISTIC_LENGTH];
        rDamageThreshold = (*(mmodel_part_old.Elements().ptr_begin()))->GetProperties()[DAMAGE_THRESHOLD];
    }
  
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void NodalDisplacementsMapping(const double& X_max,const double& X_min,const double& Y_max,const double& Y_min)
    {
        double AverageElementLength = 0.0;
        for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
        {
            AverageElementLength += (*it)->GetGeometry().Length();
        }
        AverageElementLength = AverageElementLength/mmodel_part_old.NumberOfElements();
    
        int NRows = int((Y_max-Y_min)/AverageElementLength);
        int NColumns = int((X_max-X_min)/AverageElementLength);
    
        double RowSize = (Y_max-Y_min)/NRows;
        double ColumnSize = (X_max-X_min)/NColumns;

        ElementOldCell** ElementOldMatrix;
        ElementOldMatrix = new ElementOldCell*[NRows];
        for(int i = 0; i < NRows; i++)  ElementOldMatrix[i] = new ElementOldCell[NColumns];
                
        double X_left, X_right, Y_top, Y_bot, X_me, Y_me;
        int Column_left, Column_right, Row_top, Row_bot;
        
        // Locate Old Elements in Cells
        for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
        {
            X_left = (*it)->GetGeometry().GetPoint(0).X0();
            X_right = X_left;
            Y_top = (*it)->GetGeometry().GetPoint(0).Y0();
            Y_bot = Y_top;
            for(unsigned int i = 1; i < (*it)->GetGeometry().PointsNumber(); i++)
            {
                X_me = (*it)->GetGeometry().GetPoint(i).X0();
                Y_me = (*it)->GetGeometry().GetPoint(i).Y0();
                
                if( X_me > X_right )  X_right = X_me;
                else if( X_me < X_left )  X_left = X_me;
              
                if( Y_me > Y_top ){Y_top = Y_me;}
                else if( Y_me < Y_bot ){Y_bot = Y_me;}
            }
            
            Column_left = int((X_left - X_min) / ColumnSize);
            Column_right = int((X_right - X_min) / ColumnSize);
            Row_top = int((Y_max - Y_top) / RowSize);
            Row_bot = int((Y_max - Y_bot) / RowSize);

            if(Column_left == NColumns)  Column_left = NColumns - 1;
            if(Column_right == NColumns)  Column_right = NColumns - 1;
            if(Row_top == NRows)  Row_top = NRows - 1;
            if(Row_bot == NRows)  Row_bot = NRows - 1;

            for(int i = Row_top; i <= Row_bot; i++)
            {
                for(int j = Column_left; j<= Column_right; j++)
                {
                    ElementOldMatrix[i][j].ElementOldVector.push_back((*it));
                }
            }
        }
    
        int Row, Column;
        std::vector<NodeNew*> NodeNewVector;
        
        //Locate New Nodes in Cells
        for(ModelPart::NodeIterator i = mmodel_part_new.NodesBegin(); i != mmodel_part_new.NodesEnd(); ++i)
        {
            X_me = i->X0();
            Y_me = i->Y0();

            Row = int((Y_max - Y_me) / RowSize);
            Column = int((X_me - X_min) / ColumnSize);

            if(Column == NColumns) Column = NColumns - 1;
            if(Row == NRows) Row = NRows - 1;
            
            NodeNewVector.push_back(new NodeNew((*i), Row, Column));
        }

        Element::Pointer pElementOld;
        Element::GeometryType::CoordinatesArrayType NodeLocalCoordinates;
        Vector ElementDisplacements, ElementShapeFunctions;
        double Tolerance = 1e-4;
        bool IsInside;
        
        //Locate new nodes inside old elements and interpolate displacements.
        //Triangles2D3N
        if((*mmodel_part_old.Elements().ptr_begin())->GetGeometry().PointsNumber() == 3)
        {
            ElementDisplacements = ZeroVector(3);
            for(unsigned int i = 0; i < NodeNewVector.size(); i++)
            {
                Element::GeometryType::PointType& NodeNew_i = NodeNewVector[i]->rNode;
                Row = NodeNewVector[i] -> Row;
                Column = NodeNewVector[i] -> Column;
                const Element::GeometryType::CoordinatesArrayType NodeGlobalCoordinates = NodeNew_i.Coordinates(); //Coordinates of new nodes are still in the original position
                IsInside = false;
                
                for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
                {
                    pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                    IsInside = pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates,NodeLocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                    if(IsInside)  break;                                                                        
                }
        
                if(IsInside == false) //TODO: cal??
                {
                    for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
                    {
                        pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                        pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates,NodeLocalCoordinates);
                        if(((NodeLocalCoordinates[0]+Tolerance)>=0)&&((NodeLocalCoordinates[1]+Tolerance)>=0)&&((NodeLocalCoordinates[1]-Tolerance)<=(1-NodeLocalCoordinates[0])))  break;
                    }
                }
        
                ElementShapeFunctions = this->TriangleShapeFunctions(NodeLocalCoordinates[0],NodeLocalCoordinates[1]);
        
                if( (NodeNew_i.pGetDof(DISPLACEMENT_X))->IsFixed() == false )
                {
                    ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = inner_prod(ElementShapeFunctions,ElementDisplacements);
                }
                else if( mImposedDisplacement == "Linearly_Incremented" )
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = mmodel_part_old.GetProcessInfo()[TIME_STEPS] * NodeNew_i.FastGetSolutionStepValue(IMPOSED_DISPLACEMENT)[0];
                
                if( (NodeNew_i.pGetDof(DISPLACEMENT_Y))->IsFixed() == false )
                {
                    ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[1] = inner_prod(ElementShapeFunctions,ElementDisplacements);
                }
                else if( mImposedDisplacement == "Linearly_Incremented" )
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[1] = mmodel_part_old.GetProcessInfo()[TIME_STEPS] * NodeNew_i.FastGetSolutionStepValue(IMPOSED_DISPLACEMENT)[1];
            }
        }
        else //Quadrilateral2D4N
        {
            ElementDisplacements = ZeroVector(4);
            for(unsigned int i = 0; i < NodeNewVector.size(); i++)
            {
                Element::GeometryType::PointType& NodeNew_i = NodeNewVector[i]->rNode;
                Row = NodeNewVector[i]->Row;
                Column = NodeNewVector[i]->Column;
                const Element::GeometryType::CoordinatesArrayType NodeGlobalCoordinates = NodeNew_i.Coordinates();
                IsInside = false;
        
                for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
                {
                    pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                    IsInside = pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates,NodeLocalCoordinates);
                    if(IsInside) break;
                }

                if(IsInside==false) //TODO: cal??
                {
                    for(unsigned int j = 0; j < ElementOldMatrix[Row][Column].ElementOldVector.size(); j++)
                    {
                        pElementOld = ElementOldMatrix[Row][Column].ElementOldVector[j];
                        pElementOld->GetGeometry().IsInside(NodeGlobalCoordinates,NodeLocalCoordinates);
                        if(((NodeLocalCoordinates[0]+Tolerance)>=-1)&&((NodeLocalCoordinates[0]-Tolerance)<=1)&&((NodeLocalCoordinates[1]+Tolerance)>=-1)&&((NodeLocalCoordinates[1]-Tolerance)<=1))  break;
                    }
                }
        
                ElementShapeFunctions = this->QuadrilateralShapeFunctions(NodeLocalCoordinates[0],NodeLocalCoordinates[1]);
        
                if( (NodeNew_i.pGetDof(DISPLACEMENT_X))->IsFixed() == false )
                {
                    ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    ElementDisplacements[3] = pElementOld->GetGeometry().GetPoint(3).FastGetSolutionStepValue(DISPLACEMENT)[0];
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = inner_prod(ElementShapeFunctions,ElementDisplacements);
                }
                else if( mImposedDisplacement == "Linearly_Incremented" )
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[0] = mmodel_part_old.GetProcessInfo()[TIME_STEPS]*NodeNew_i.FastGetSolutionStepValue(IMPOSED_DISPLACEMENT)[0];

                if( (NodeNew_i.pGetDof(DISPLACEMENT_Y))->IsFixed()==false )
                {
                    ElementDisplacements[0] = pElementOld->GetGeometry().GetPoint(0).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    ElementDisplacements[1] = pElementOld->GetGeometry().GetPoint(1).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    ElementDisplacements[2] = pElementOld->GetGeometry().GetPoint(2).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    ElementDisplacements[3] = pElementOld->GetGeometry().GetPoint(3).FastGetSolutionStepValue(DISPLACEMENT)[1];
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[1] = inner_prod(ElementShapeFunctions,ElementDisplacements);
                }
                else if( mImposedDisplacement == "Linearly_Incremented" )
                    NodeNew_i.FastGetSolutionStepValue(DISPLACEMENT)[1] = mmodel_part_old.GetProcessInfo()[TIME_STEPS]*NodeNew_i.FastGetSolutionStepValue(IMPOSED_DISPLACEMENT)[1];
            }
        }
    
        //Deallocate memory
        for(int i=0;i<NRows;i++)
            delete[] ElementOldMatrix[i];
        delete[] ElementOldMatrix;
    
        for(unsigned int i=0;i<NodeNewVector.size(); i++)
            delete NodeNewVector[i];
    
        std::cout << "Nodal Displacements Mapped" << std::endl;
    
    }//NodalDisplacementsMapping_End

    //------------------------------------------------------------------------------------

    Vector TriangleShapeFunctions(const double& rx_local, const double& ry_local)
    {
        Vector TSF = ZeroVector(3);
        TSF[0] = 1-rx_local-ry_local;
        TSF[1] = rx_local;
        TSF[2] = ry_local;
        return TSF;    
    }

    //------------------------------------------------------------------------------------
    
    Vector QuadrilateralShapeFunctions(const double& rx_local, const double& ry_local)
    {
        Vector QSF = ZeroVector(4);
        QSF[0] = (1-rx_local-ry_local+rx_local*ry_local)/4;
        QSF[1] = (1+rx_local-ry_local-rx_local*ry_local)/4;
        QSF[2] = (1+rx_local+ry_local+rx_local*ry_local)/4;
        QSF[3] = (1-rx_local+ry_local-rx_local*ry_local)/4;
        return QSF;    
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GaussPointStateVariableMapping(const double& X_max,const double& X_min,const double& Y_max,const double& Y_min,const double& rCharacteristicLength,const double& rDamageThreshold)
    {
        std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector;
        int Row, Column;
        double X_me, Y_me, X_other, Y_other, Distance;
        ConstitutiveLaw::Pointer Me;
        ConstitutiveLaw::Pointer Other;
        Vector Trian_GPLocalCoord = ZeroVector(3);
        std::vector<Vector> Quad_GPLocalCoordVector(4);
        Element::GeometryType::CoordinatesArrayType GPGlobalCoord;
        std::vector<GaussPointNew*> GaussPointNewVector;

        int NRows = int((Y_max - Y_min) / (rCharacteristicLength * 2));
        int NColumns = int((X_max - X_min) / (rCharacteristicLength * 2));
    
        double RowSize = (Y_max - Y_min) / NRows;
        double ColumnSize = (X_max - X_min) / NColumns;

        GaussPointOldCell** pGaussPointOldMatrix;
        pGaussPointOldMatrix = new GaussPointOldCell*[NRows];
        for(int i=0;i<NRows;i++)  pGaussPointOldMatrix[i] = new GaussPointOldCell[NColumns];
        
        //Locate old Gauss points in cells
        //Triangles2D3N
        if((*mmodel_part_old.Elements().ptr_begin())->GetGeometry().PointsNumber() == 3)
        {
            // GP local coordinates
            Trian_GPLocalCoord[0] = 1 / 3;
            Trian_GPLocalCoord[1] = 1 / 3;
            const Element::GeometryType::CoordinatesArrayType GPLocalCoord = Trian_GPLocalCoord;
            
            for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
            {
                (*it)->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER,ConstitutiveLawVector,mmodel_part_old.GetProcessInfo());

                (*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);

                X_me = GPGlobalCoord[0];
                Y_me = GPGlobalCoord[1];

                Row    = int((Y_max - Y_me) / RowSize);
                Column = int((X_me - X_min) / ColumnSize);

                if(Row == NRows) Row = NRows - 1;
                if(Column == NColumns) Column = NColumns-1;

                pGaussPointOldMatrix[Row][Column].GaussPointOldVector.push_back(GaussPointOld(ConstitutiveLawVector[0], X_me, Y_me));
            }
        }
        else //Quadrilateral2D4N
        {
            // GP local coordinates
            Quad_GPLocalCoordVector[0] = ZeroVector(3);
            Quad_GPLocalCoordVector[0][0] = -sqrt(3)/3;
            Quad_GPLocalCoordVector[0][1] = -sqrt(3)/3;
            Quad_GPLocalCoordVector[1] = ZeroVector(3);
            Quad_GPLocalCoordVector[1][0] = sqrt(3)/3;
            Quad_GPLocalCoordVector[1][1] = -sqrt(3)/3;
            Quad_GPLocalCoordVector[2] = ZeroVector(3);
            Quad_GPLocalCoordVector[2][0] = sqrt(3)/3;
            Quad_GPLocalCoordVector[2][1] = sqrt(3)/3;
            Quad_GPLocalCoordVector[3] = ZeroVector(3);
            Quad_GPLocalCoordVector[3][0] = -sqrt(3)/3;
            Quad_GPLocalCoordVector[3][1] = sqrt(3)/3;
      
            for(ElementsArrayType::ptr_iterator it = mmodel_part_old.Elements().ptr_begin(); it != mmodel_part_old.Elements().ptr_end(); ++it)
            {
                (*it)->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER,ConstitutiveLawVector,mmodel_part_old.GetProcessInfo());

                for(unsigned int i = 0; i < 4; i++)
                {
                    const Element::GeometryType::CoordinatesArrayType GPLocalCoord = Quad_GPLocalCoordVector[i];
            
                    (*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);
            
                    X_me = GPGlobalCoord[0];
                    Y_me = GPGlobalCoord[1];
            
                    Row = int((Y_max-Y_me)/RowSize);
                    Column = int((X_me-X_min)/ColumnSize);
                    
                    if(Row==NRows) Row = NRows-1;
                    if(Column==NColumns) Column = NColumns-1;
                    
                    pGaussPointOldMatrix[Row][Column].GaussPointOldVector.push_back(GaussPointOld(ConstitutiveLawVector[i],X_me,Y_me));
                }
            }
        }
    
        //Locate new Gauss points in cells
        //Triangles2D3N
        if((*mmodel_part_new.Elements().ptr_begin())->GetGeometry().PointsNumber()==3)
        {
            const Element::GeometryType::CoordinatesArrayType GPLocalCoord = Trian_GPLocalCoord;
            
            for(ElementsArrayType::ptr_iterator it = mmodel_part_new.Elements().ptr_begin(); it != mmodel_part_new.Elements().ptr_end(); ++it)
            {
                (*it)->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER,ConstitutiveLawVector,mmodel_part_new.GetProcessInfo());
                
                (*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);
                
                X_me = GPGlobalCoord[0];
                Y_me = GPGlobalCoord[1];

                Row = int((Y_max - Y_me) / RowSize);
                Column = int((X_me - X_min) / ColumnSize);

                if(Row == NRows) Row = NRows - 1;
                if(Column == NColumns) Column = NColumns - 1;

                GaussPointNewVector.push_back(new GaussPointNew(ConstitutiveLawVector[0],X_me, Y_me,Row,Column));
            }
        }
        else //Quadrilateral2D4N
        {
            for(ElementsArrayType::ptr_iterator it = mmodel_part_new.Elements().ptr_begin(); it != mmodel_part_new.Elements().ptr_end(); ++it)
            {
                (*it)->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW_POINTER,ConstitutiveLawVector,mmodel_part_new.GetProcessInfo());

                for(unsigned int i = 0; i < 4; i++)
                { 
                    const Element::GeometryType::CoordinatesArrayType GPLocalCoord = Quad_GPLocalCoordVector[i];
            
                    (*it)->GetGeometry().GlobalCoordinates(GPGlobalCoord,GPLocalCoord);

                    X_me = GPGlobalCoord[0];
                    Y_me = GPGlobalCoord[1];

                    Row = int((Y_max-Y_me)/RowSize);
                    Column = int((X_me-X_min)/ColumnSize);

                    if(Row==NRows) Row = NRows-1;
                    if(Column==NColumns) Column = NColumns-1;

                    GaussPointNewVector.push_back(new GaussPointNew(ConstitutiveLawVector[i],X_me, Y_me,Row,Column));
                }
            }
        }   
    
        //Transfer state variables from old Gauss points to new Gauss Points (nonlocal average)
        double IntegrationCoefficient, StateVariable, Numerator, WeightingFunctionDenominator;
        for(unsigned int i=0; i < GaussPointNewVector.size(); i++)
        {
            Me = GaussPointNewVector[i]->pConstitutiveLaw;
            X_me = GaussPointNewVector[i]->X_coord;
            Y_me = GaussPointNewVector[i]->Y_coord;
            Row = GaussPointNewVector[i]->Row;
            Column = GaussPointNewVector[i]->Column;
            Numerator = 0.0;
            WeightingFunctionDenominator = 0.0;
            
            //Search in my cell
            for(unsigned int j=0; j<pGaussPointOldMatrix[Row][Column].GaussPointOldVector.size(); j++)
            {
                Other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].pConstitutiveLaw;
                X_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].X_coord;
                Y_other = pGaussPointOldMatrix[Row][Column].GaussPointOldVector[j].Y_coord;
        
                Distance = sqrt((X_other-X_me)*(X_other-X_me) + (Y_other-Y_me)*(Y_other-Y_me));

                if(Distance <= rCharacteristicLength) //TODO: es podria calcular amb tots els de la teva cel·la...
                {
                    IntegrationCoefficient = Other->GetValue(INTEGRATION_COEFFICIENT,IntegrationCoefficient);
                    StateVariable = Other->GetValue(STATE_VARIABLE,StateVariable);

                    Numerator += IntegrationCoefficient*exp(-4*Distance*Distance/(rCharacteristicLength*rCharacteristicLength))*StateVariable;
                    WeightingFunctionDenominator += IntegrationCoefficient*exp(-4*Distance*Distance/(rCharacteristicLength*rCharacteristicLength));
                }
            }
            
            //Search in adjacent cells
            if(sqrt((X_min+ColumnSize*(Column+1)-X_me)*(X_min+ColumnSize*(Column+1)-X_me) + (Y_max-RowSize*(Row+1)-Y_me)*(Y_max-RowSize*(Row+1)-Y_me)) < rCharacteristicLength)
            {
                if((Row+1)<NRows)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row+1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if((Column+1)<NColumns)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column+1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if(((Row+1)<NRows) && ((Column+1)<NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row+1][Column+1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
            }
            else if( sqrt((X_min+ColumnSize*(Column)-X_me)*(X_min+ColumnSize*(Column)-X_me) + (Y_max-RowSize*(Row+1)-Y_me)*(Y_max-RowSize*(Row+1)-Y_me)) < rCharacteristicLength)
            {
                if((Row+1)<NRows)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row+1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if((Column-1)>=0)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column-1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if(((Row+1)<NRows) && ((Column-1)>=0))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row+1][Column-1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
            }
            else if( sqrt((X_min+ColumnSize*(Column)-X_me)*(X_min+ColumnSize*(Column)-X_me) + (Y_max-RowSize*(Row)-Y_me)*(Y_max-RowSize*(Row)-Y_me)) < rCharacteristicLength)
            {
                if((Row-1)>=0)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row-1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if((Column-1)>=0)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column-1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if(((Row-1)>=0) && ((Column-1)>=0))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row-1][Column-1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
            }
            else if( sqrt((X_min+ColumnSize*(Column+1)-X_me)*(X_min+ColumnSize*(Column+1)-X_me) + (Y_max-RowSize*(Row)-Y_me)*(Y_max-RowSize*(Row)-Y_me)) < rCharacteristicLength)
            {
                if((Row-1)>=0)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row-1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if((Column+1)<NColumns)  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column+1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                if(((Row-1)>=0) && ((Column+1)<NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row-1][Column+1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
            }
            else
            {
                if((int((X_me-X_min+rCharacteristicLength)/ColumnSize)>Column) && ((Column+1)<NColumns))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column+1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                else if((int((X_me-X_min-rCharacteristicLength)/ColumnSize)<Column) && ((Column-1)>=0))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row][Column-1],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
        
                if((int((Y_max-Y_me+rCharacteristicLength)/RowSize)>Row) && ((Row+1)<NRows))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row+1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
                else if((int((Y_max-Y_me-rCharacteristicLength)/RowSize)<Row) && ((Row-1)>=0))  this->SearchInAdjacentCell(GaussPointNewVector[i],pGaussPointOldMatrix[Row-1][Column],Numerator,WeightingFunctionDenominator,rCharacteristicLength);
            }
      
            if(fabs(WeightingFunctionDenominator)<1e-15)
                StateVariable = rDamageThreshold;
            else  
                StateVariable = Numerator/WeightingFunctionDenominator;
      
            Me->SetValue(STATE_VARIABLE,StateVariable,mmodel_part_new.GetProcessInfo());
            Me->SetValue(STATE_VARIABLE_EQUILIBRIUM,StateVariable,mmodel_part_new.GetProcessInfo());
        }

        //Deallocate memory
        for(int i=0;i<NRows;i++)
            delete [] pGaussPointOldMatrix[i];
        delete [] pGaussPointOldMatrix;
    
        for(unsigned int i=0;i<GaussPointNewVector.size(); i++) 
            delete GaussPointNewVector[i];
    
        std::cout << "Gauss Point State Variable Mapped" << std::endl;
    
    }//GaussPointStateVariableMapping_End

    //------------------------------------------------------------------------------------
    
    void SearchInAdjacentCell(const GaussPointNew* pGaussPointNew, const GaussPointOldCell& NeighbourCell, double& rNumerator, double& rWeightingFunctionDenominator, const double& rCharacteristicLength)
    {
        double X_me = pGaussPointNew->X_coord;
        double Y_me = pGaussPointNew->Y_coord;
    
        ConstitutiveLaw::Pointer Other;
        double X_other,Y_other,Distance;
        double IntegrationCoefficient, StateVariable;
    
        for(unsigned int j=0; j<NeighbourCell.GaussPointOldVector.size(); j++)
        {
            Other = NeighbourCell.GaussPointOldVector[j].pConstitutiveLaw;
            X_other = NeighbourCell.GaussPointOldVector[j].X_coord;
            Y_other = NeighbourCell.GaussPointOldVector[j].Y_coord;

            Distance = sqrt((X_other-X_me)*(X_other-X_me) + (Y_other-Y_me)*(Y_other-Y_me));
      
            if(Distance <= rCharacteristicLength) //TODO: es podria calcular amb tots els de les cel·les vehines...
            {
                IntegrationCoefficient = Other->GetValue(INTEGRATION_COEFFICIENT,IntegrationCoefficient);
                StateVariable = Other->GetValue(STATE_VARIABLE,StateVariable);

                rNumerator += IntegrationCoefficient*exp(-4*Distance*Distance/(rCharacteristicLength*rCharacteristicLength))*StateVariable;
                rWeightingFunctionDenominator += IntegrationCoefficient*exp(-4*Distance*Distance/(rCharacteristicLength*rCharacteristicLength));
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void TransferProcessInfoVariables()
    {
        //Arc-length parameters
        mmodel_part_new.GetProcessInfo()[LOAD_FACTOR] = mmodel_part_old.GetProcessInfo()[LOAD_FACTOR];
        mmodel_part_new.GetProcessInfo()[RADIUS_FACTOR] = mmodel_part_old.GetProcessInfo()[RADIUS_FACTOR];
        
        //To compute linearly incremented loads
        mmodel_part_new.GetProcessInfo()[TIME_STEPS] = mmodel_part_old.GetProcessInfo()[TIME_STEPS];
    }

};//Class

} /* namespace Kratos.*/

#endif /* KRATOS_MAPPING_VARIABLES_PROCESS defined */
