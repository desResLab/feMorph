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
    std::vector<int> elementConnections;
    std::vector<int> elementFaces;
    // Constructor and Destructor
    femElement(int number, int prop, int totalNodes, int* connections);
    femElement(const femElement* other);
    virtual ~femElement();
    // Member Functions
    void evalElementCentroid(std::vector<femNode*> &nodeList, double* centroid);
    double evalPointToElementDistance(double* pointCoords, std::vector<femNode*> &nodeList);
    void InterpolateElementDisplacements(double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps);
    virtual bool isNodeInsideElement(double* pointCoords,std::vector<femNode*> &nodeList){return false;}
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){}
    void CheckRandomWalkingCriterion(int localFaceID, double* nodeCoords,std::vector<femFace*> &faceList,std::vector<femNode*> &nodeList,bool &isOnFace, bool &isOnOppositeSide);
    int getAdjacentElement(int localFaceID, std::vector<femFace*> &faceList);
};

// TETRAHEDRAL ELEMENT
class femTetra4: public femElement
{
  public:
    femTetra4(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){}
    femTetra4(const femElement* other):femElement(other){}
    virtual ~femTetra4(){}
    // Member Functions
    bool isNodeInsideElement(double* pointCoords,std::vector<femNode*> &nodeList);
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    void AssembleTetCoordsMat(std::vector<femNode*> &nodeList, double** coordMat);

};

// TETRAHEDRAL ELEMENT
class femTetra10: public femElement
{
  public:
    femTetra10(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){}
    femTetra10(const femElement* other):femElement(other){}
    virtual ~femTetra10(){}
    // Member Functions
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    virtual bool isNodeInsideElement(double* pointCoords,std::vector<femNode*> &nodeList);
};

// TETRAHEDRAL ELEMENT
class femHexa8: public femElement
{
  public:
    femHexa8();
    virtual ~femHexa8(){}
    // Member Functions
    virtual bool isNodeInsideElement(double* nodeCoords);
};

#endif // FEMELEMENT_H
