# include "math.h"
# include "femRod.h"
# include "femException.h"

femRod::femRod(int number, int prop, int totalNodes, int* connections, double currArea):femElement(number,prop,totalNodes,connections){
  numberOfNodes = 2;
  area = currArea;
  dims = d1;
}

// Member Functions
bool femRod::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

// Eval Rod Volume
double femRod::EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){
  // Eval Rod Length
  int node1 = elementConnections[0];
  int node2 = elementConnections[1];
  double dx = nodeList[node2]->coords[0] - nodeList[node1]->coords[0];
  double dy = nodeList[node2]->coords[1] - nodeList[node1]->coords[1];
  double dz = nodeList[node2]->coords[2] - nodeList[node1]->coords[2];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

// Positive Volume Evaluation
double femRod::EvalMixProduct(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

// Fix Connectivities by Swapping End Nodes
void femRod::fixConnectivities(std::vector<femNode*> &nodeList){
  int temp = 0;
  temp = elementConnections[0];
  elementConnections[0] = elementConnections[1];
  elementConnections[1] = temp;
}

// Eval Local shape functions
void femRod::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  // Insert the two shape functions
  shapeFunction.clear();
  shapeFunction.push_back(0.5*(1.0 - coord1));
  shapeFunction.push_back(0.5*(1.0 + coord1));
}

// Eval Local shape functions derivatives
void femRod::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  // Insert the two shape functions
  femDoubleVec temp;
  shapeDeriv.clear();
  temp.push_back(-0.5);
  temp.push_back(0.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
  temp.clear();
  temp.push_back(0.5);
  temp.push_back(0.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
}

