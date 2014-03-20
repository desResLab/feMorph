#ifndef FEMELEMENT_H
#define FEMELEMENT_H

#include <vector>

#include "femNode.h"
#include "femFace.h"

// GENERIC ELEMENT
class femElement
{
  public:
    // Data Members
    int elementNumber;
    int propertyNumber;
    int numberOfNodes;
    std::vector<int> elementConnections;
    std::vector<int> elementFaces;
    std::vector<int> elementEdges;
    // Constructor and Destructor
    femElement(int number, int prop, int totalNodes, int* connections);
    femElement(const femElement* other);
    virtual ~femElement();
    // Member Functions
    // Virtual
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){return false;}
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){}
    virtual bool is2D(){return false;}
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){return 0.0;}
    // Not Virtual
    void   evalElementCentroid(std::vector<femNode*> &nodeList, double* centroid);
    double evalPointToElementDistance(double* pointCoords, std::vector<femNode*> &nodeList);
    void   InterpolateElementDisplacements(double dispFactor, double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps);
    void   CheckRandomWalkingCriterion(int localFaceID, double* nodeCoords,std::vector<femFace*> &faceList,std::vector<femNode*> &nodeList,bool &isOnFace, bool &isOnOppositeSide);
    int    getAdjacentElement(int localFaceID, std::vector<femFace*> &faceList);
    void   CreateBoundingBoxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &boxNodeList);
    void   CreateMinMaxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &minMaxNodeList);
};

// TETRAHEDRAL ELEMENT
class femTetra4: public femElement
{
  public:
    femTetra4(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 4;}
    femTetra4(const femElement* other):femElement(other){numberOfNodes = 4;}
    virtual ~femTetra4(){}
    // Member Functions
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual bool is2D(){return false;}
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList);
    virtual double EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList);
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    void AssembleTetCoordsMat(double dispFactor, std::vector<femNode*> &nodeList, double** coordMat);
};

// TETRAHEDRAL ELEMENT
class femTetra10: public femElement
{
  public:
    femTetra10(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 10;}
    femTetra10(const femElement* other):femElement(other){numberOfNodes = 10;}
    virtual ~femTetra10(){}
    // Member Functions
    virtual bool is2D(){return false;}
    virtual void EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual bool isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual double EvalVolume(std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList);
};

// TETRAHEDRAL ELEMENT
class femHexa8: public femElement
{
  public:
    femHexa8(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 8;}
    virtual ~femHexa8(){}
    // Member Functions
    virtual bool is2D(){return false;}
    virtual bool isNodeInsideElement(double* nodeCoords);
    virtual double EvalVolume(std::vector<femNode*> &nodeList){return 0.0;}
};

// TRIANGULAR ELEMENT
class femTri3: public femElement
{
  public:
    femTri3(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){numberOfNodes = 3;}
    virtual ~femTri3(){}
    // Member Functions
    virtual bool is2D(){return true;}
    virtual double EvalVolume(std::vector<femNode*> &nodeList){return 0.0;}
    virtual double EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList);
};


#endif // FEMELEMENT_H
