# include "femElement.h"
# include "femException.h"
# include "femConstants.h"
# include "femUtils.h"

// EVAL SHAPE FUNCTION AT GIVEN LOCATION
void femHexa8::evalShapeFunction(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
 shapeFunction.clear();
 shapeFunction.reserve(numberOfNodes);
 double c1[8] = {-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
 double c2[8] = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
 double c3[8] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
 double currN = 0.0;
 for(int loopA=0;loopA<numberOfNodes;loopA++){
   currN = (1.0/8.0)*(1.0 + coord1 * c1[loopA])*
                     (1.0 + coord2 * c2[loopA])*
                     (1.0 + coord3 * c3[loopA]);
   shapeFunction.push_back(currN);
 }
}

void femHexa8::evalLocalShapeFunctionDerivative(std::vector<femNode*> &nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  shapeDeriv.clear();
  shapeDeriv.reserve(numberOfNodes);
  double c1[8] = {-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
  double c2[8] = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
  double c3[8] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
  femDoubleVec temp;
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    temp.clear();
    temp.push_back((1.0/8.0)*c1[loopA]*(1.0 + coord2 * c2[loopA])*(1.0 + coord3 * c3[loopA]));
    temp.push_back((1.0/8.0)*c2[loopA]*(1.0 + coord1 * c1[loopA])*(1.0 + coord3 * c3[loopA]));
    temp.push_back((1.0/8.0)*c3[loopA]*(1.0 + coord1 * c1[loopA])*(1.0 + coord2 * c2[loopA]));
    shapeDeriv.push_back(temp);
  }
}

bool femHexa8::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

double femHexa8::EvalVolume(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

// EVAL ELEMENT MIXED PRODUCT
double femHexa8::EvalMixProduct(std::vector<femNode*> &nodeList){
  // Init
  double vecA[3] = {0.0};
  double vecB[3] = {0.0};
  double vecC[3] = {0.0};
  double vec1[3] = {0.0};
  double currProd = 0.0;
  // Get three vectors
  vecA[0] = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  vecA[1] = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  vecA[2] = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  vecB[0] = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  vecB[1] = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  vecB[2] = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  vecC[0] = nodeList[elementConnections[4]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  vecC[1] = nodeList[elementConnections[4]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  vecC[2] = nodeList[elementConnections[4]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  // Get external product
  femUtils::Do3DExternalProduct(vecA,vecB,vec1);
  // Get Internal Product
  currProd = femUtils::Do3DInternalProduct(vec1,vecC);
  return currProd;
}

// SWAP NODES
void femHexa8::fixConnectivities(std::vector<femNode*> &nodeList){
  // Swap 2 with 3
  int temp = elementConnections[2];
  elementConnections[2] = elementConnections[3];
  elementConnections[3] = temp;

  // Swap 6 with 7
  temp = elementConnections[6];
  elementConnections[6] = elementConnections[7];
  elementConnections[7] = temp;
}



