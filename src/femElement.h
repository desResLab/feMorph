#ifndef FEMELEMENT_H
#define FEMELEMENT_H

# include <vector>

# include "femNode.h"
# include "femFace.h"
# include "femTypes.h"
# include "femIntegrationRule.h"
# include "femTime.h"

// GENERIC ELEMENT
class femElement{
  public:
    // Data Members
    int elementNumber;
    int propertyNumber;
    int numberOfNodes;
    elDim dims;

    std::vector<int> elementConnections;
    std::vector<int> elementFaces;
    std::vector<int> elementEdges;
    // Constructor and Destructor
    femElement(int number, int prop, int totalNodes, int* connections);
    femElement(const femElement* other);
    virtual ~femElement();
    // Member Functions
    // Virtual
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList);
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);

    // Common
    void   evalJacobianMatrix(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, elDim dims, femDoubleMat &jacMat);
    double evalJacobian(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3);
    void   evalGlobalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, double& detJ, femDoubleMat &globShDeriv);
    void   evalGeometricMatrix(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3,femDoubleMat& elGeomMat);

    // Temporary
    void getBCCoords(double bcCoord1, double bcCoord2, double bcCoord3,
                     femIntVec conn,
                     double& newCoord1, double& newCoord2, double& newCoord3);


    // Not Virtual
    void   evalElementCentroid(std::vector<femNode*> &nodeList, double* centroid);
    double evalPointToElementDistance(double* pointCoords, std::vector<femNode*> &nodeList);
    void   InterpolateElementDisplacements(double dispFactor, double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps);
    void   CheckRandomWalkingCriterion(int localFaceID, double* nodeCoords,std::vector<femFace*> &faceList,std::vector<femNode*> &nodeList,bool &isOnFace, bool &isOnOppositeSide);
    int    getAdjacentElement(int localFaceID, std::vector<femFace*> &faceList);
    void   CreateBoundingBoxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &boxNodeList);
    void   CreateMinMaxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &minMaxNodeList);    
    double checkMinDetJ(std::vector<femNode*> &nodeList, femIntegrationRule* rule);

    // BOUNDARY CONDITIONS
    femElement* createBoudaryFace(int faceID);

    // POISSON PROBLEM
    void formPoissonMatrix(std::vector<femNode*>& nodeList,femIntegrationRule* rule,femDoubleVec diffusivity,femDoubleMat &elMat);
    void formPoissonSource(std::vector<femNode*>& nodeList,femIntegrationRule* rule, double sourceValue,femDoubleVec &elSourceVec);
    void formPoissonNeumannBC();

    // INCOMPRESSIBLE NAVIER STOKES
    void formNS_LHS(std::vector<femNode*> nodeList,
                    femIntegrationRule* rule,
                    double rho,
                    double viscosity,
                    femDoubleMat solution,
                    int schemeType,
                    int nodeDOFs,
                    double timeStep,
                    double alphaM,
                    double alphaF,
                    double gamma,
                    femDoubleBlockMat &elMat);

    void formNS_RHS(std::vector<femNode*> &nodeList,
                    femIntegrationRule* rule,
                    double sourceValue,
                    femDoubleVec diffusivity,
                    femDoubleMat solution,
                    femDoubleMat solution_Dot,
                    double timeStep,
                    double alphaM, double alphaF, double gamma,
                    femDoubleBlockVec &elRhs);

    // ADVECTION DIFFUSION
    // Steady State
    void formAdvDiffLHS(std::vector<femNode*> nodeList,femIntegrationRule* rule,femDoubleVec diffusivity, femDoubleVec velocity,int schemeType,femDoubleMat &elMat);
    void formAdvDiffRHS(std::vector<femNode*> nodeList,femIntegrationRule* rule,double sourceValue,femDoubleVec diffusivity,femDoubleVec velocity,int schemeType,femDoubleVec &elRhs);
    void formWeakBC(std::vector<femNode*> nodeList,
                    femElement* parentElement,
                    femIntegrationRule* rule,
                    femDoubleVec diffusivity,femDoubleVec velocity,femDoubleVec elNormal, double elBCValue,
                    femDoubleMat &elMat,femDoubleVec &elVec);
    // Transient
    void formTransientAdvDiffLHS(std::vector<femNode*> &nodeList,
                                 femIntegrationRule* rule,
                                 femDoubleVec diffusivity,
                                 femDoubleMat solution,
                                 double timeStep,
                                 double alphaM, double alphaF, double gamma,
                                 femDoubleMat &elMat);

    void formTransientAdvDiffRHS(std::vector<femNode*> &nodeList,
                                 femIntegrationRule* rule,
                                 double sourceValue,
                                 femDoubleVec diffusivity,
                                 int variableID,
                                 femDoubleMat solution,
                                 femDoubleMat solution_Dot,
                                 double timeStep,
                                 double alphaM, double alphaF, double gamma,
                                 femDoubleVec &elRhs);


    void assembleMass(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double** massMat);
    void assembleStiffness(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double diffusivity, femDoubleMat &stiffnessMat);


    // INTEGRATE NODAL VECTOR
    double integrateNodalVector(std::vector<femNode*> nodeList,femIntegrationRule* rule,femDoubleVec nodeVec);
};

// TETRAHEDRAL ELEMENT
class femTetra4: public femElement{
  public:
    femTetra4(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 4;dims = d3;}
    femTetra4(const femElement* other):femElement(other){numberOfNodes = 4;}
    virtual ~femTetra4(){}
    // Member Functions
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList);
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual bool   isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);
    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);

    void AssembleTetCoordsMat(double dispFactor, std::vector<femNode*> &nodeList, double** coordMat);
};

// TETRAHEDRAL ELEMENT
class femTetra10: public femElement{
  public:
    femTetra10(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 10;dims = d3;}
    femTetra10(const femElement* other):femElement(other){numberOfNodes = 10;}
    virtual ~femTetra10(){}

    // Member Functions
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);
};

// TETRAHEDRAL ELEMENT
class femHexa8: public femElement{
  public:
    femHexa8(int number, int prop, int totalNodes, int* connections);
    virtual ~femHexa8(){}

    // Member Functions
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList);
    // Positive Volume Evaluation
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);
};

// TRIANGULAR ELEMENT
class femTri3: public femElement{
  public:
    femTri3(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 3;dims = d2;}
    virtual ~femTri3(){}

    // Member Functions
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);
};

// QUAD4 ELEMENT
class femQuad4: public femElement{
  public:
    femQuad4(int number, int prop, int totalNodes, int* connections);
    virtual ~femQuad4(){}

    // Member Functions
    virtual bool   is2D(){return true;}
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList){return 0.0;}
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);
};

#endif // FEMELEMENT_H
