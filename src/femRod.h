#ifndef FEMROD_H
#define FEMROD_H

# include "femElement.h"

// ROD ELEMENT
class femRod: public femElement{
  public:
    // Data
    double area;
    // Constructor/Descructor
    femRod(int number, int prop, int totalNodes, int* connections,double currArea);
    virtual ~femRod(){}

    // Member Functions
    virtual bool   isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList);
    virtual double EvalVolume(double dispFactor, std::vector<femNode*> &nodeList);
    // Positive Volume Evaluation
    virtual double EvalMixProduct(std::vector<femNode*> &nodeList);
    virtual void   fixConnectivities(std::vector<femNode*> &nodeList);

    // Element Calculation
    virtual void   evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction);
    virtual void   evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv);
};

#endif // FEMROD_H
