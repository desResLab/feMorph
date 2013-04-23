#ifndef FEMELEMENT_H
#define FEMELEMENT_H

#include <vector>

#include "femNode.h"

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
    femElement(femElement* other);
    virtual ~femElement();
    // Member Functions
    void evalCentroid(double* centroid);
    virtual bool isNodeInsideElement(double* nodeCoords){return false;}
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){}
    virtual void InterpolateElementDisplacements(double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps);

};

// TETRAHEDRAL ELEMENT
class femTetra4: public femElement
{
  public:
    femTetra4(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){}
    virtual ~femTetra4(){}
    // Member Functions
    virtual bool isNodeInsideElement(double* pointCoords,std::vector<femNode*> &nodeList);
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
    void AssembleTetCoordsMat(std::vector<femNode*> &nodeList, double** coordMat);

};

// TETRAHEDRAL ELEMENT
class femTetra10: public femElement
{
  public:
    femTetra10();
    virtual ~femTetra10(){};
    // Member Functions
    virtual bool isNodeInsideElement(double* nodeCoords);
    virtual void EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords);
};

// TETRAHEDRAL ELEMENT
class femHexa8: public femElement
{
  public:
    femHexa8();
    virtual ~femHexa8(){};
    // Member Functions
    virtual bool isNodeInsideElement(double* nodeCoords);
};

#endif // FEMELEMENT_H
