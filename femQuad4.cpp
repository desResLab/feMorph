# include "femElement.h"
# include "femException.h"

// CONSTRUCTOR
femQuad4::femQuad4(int number, int prop, int totalNodes, int* connections):femElement(number,prop,totalNodes,connections){
  numberOfNodes = 4;
  dims = d2;
}

void femQuad4::fixConnectivities(std::vector<femNode*> &nodeList){
  // Invert local Node 1 and 2
  int temp = elementConnections[1];
  elementConnections[1] = elementConnections[2];
  elementConnections[2] = temp;
}
void femQuad4::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  shapeFunction.clear();
  shapeFunction.reserve(numberOfNodes);
  double c1[4] = {-1.0,-1.0,+1.0,+1.0};
  double c2[4] = {-1.0,+1.0,+1.0,-1.0};
  double currN = 0.0;
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    currN = (1.0/4.0)*(1.0 + coord1 * c1[loopA])*
                      (1.0 + coord2 * c2[loopA]);
    shapeFunction.push_back(currN);
  }
}
void femQuad4::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  shapeDeriv.clear();
  shapeDeriv.reserve(numberOfNodes);
  double c1[4] = {-1.0,-1.0,+1.0,+1.0};
  double c2[4] = {-1.0,+1.0,+1.0,-1.0};
  femDoubleVec temp;
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    temp.clear();
    temp.push_back((1.0/4.0)*c1[loopA]*(1.0 + coord2 * c2[loopA]));
    temp.push_back((1.0/4.0)*c2[loopA]*(1.0 + coord1 * c1[loopA]));
    temp.push_back(0.0);
    shapeDeriv.push_back(temp);
  }
}

