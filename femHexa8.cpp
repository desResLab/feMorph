# include "femElement.h"
# include "femException.h"
# include "femConstants.h"
# include "femUtils.h"

// EVAL JACOBIAN MATRIX
void evalJacobianMatrixLocal(int numberOfNodes, femDoubleMat elNodeCoords, femDoubleMat shLocalDerivs, femDoubleMat &jacMat){
  // Assemble Jacobian Matrix
  for(int loopA=0;loopA<kDims;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      jacMat[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        jacMat[loopA][loopB] += shLocalDerivs[loopC][loopA] * elNodeCoords[loopC][loopB];
      }
    }
  }
}

void evalShLocalDerivative(double coord1, double coord2, double coord3, femDoubleMat &shLocalDerivs){
  shLocalDerivs.clear();
  shLocalDerivs.reserve(8);
  double c1[8] = {-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
  double c2[8] = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
  double c3[8] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
  femDoubleVec temp;
  for(int loopA=0;loopA<8;loopA++){
    temp.clear();
    temp.push_back((1.0/8.0)*c1[loopA]*(1.0 + coord2 * c2[loopA])*(1.0 + coord3 * c3[loopA]));
    temp.push_back((1.0/8.0)*c2[loopA]*(1.0 + coord1 * c1[loopA])*(1.0 + coord3 * c3[loopA]));
    temp.push_back((1.0/8.0)*c3[loopA]*(1.0 + coord1 * c1[loopA])*(1.0 + coord2 * c2[loopA]));
    shLocalDerivs.push_back(temp);
  }
}

bool femHexa8::is2D(){
  return false;
}

bool femHexa8::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

double femHexa8::EvalVolume(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

void femHexa8::swapNodes(){
  throw femException("Not Implemented.\n");
}

// EVAL SHAPE FUNCTION AT GIVEN LOCATION
void femHexa8::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
 shapeFunction.clear();
 shapeFunction.reserve(8);
 double c1[8] = {-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0};
 double c2[8] = {-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,1.0};
 double c3[8] = {-1.0,-1.0,-1.0,-1.0,1.0,1.0,1.0,1.0};
 double currN = 0.0;
 for(int loopA=0;loopA<8;loopA++){
   currN = (1.0/8.0)*(1.0 + coord1 * c1[loopA])*
                     (1.0 + coord2 * c2[loopA])*
                     (1.0 + coord3 * c3[loopA]);
   shapeFunction.push_back(currN);
 }
}
// EVAL SHAPE FUNCTION DERIVATIVES RESPECT TO GLOBAL COORDINATES
void femHexa8::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){

  // Compute Node Coordinate Vector
  int currNode = 0;
  femDoubleMat elNodeCoords;
  elNodeCoords.resize(8);
  for(int loopA=0;loopA<8;loopA++){
    elNodeCoords[loopA].resize(kDims);
  }
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    currNode = elementConnections[loopA];
    for(int loopB=0;loopB<kDims;loopB++){
      elNodeCoords[loopA][loopB] = nodeList[currNode]->coords[loopB];
    }
  }

  // Compute Local Derivatives
  femDoubleMat shLocalDerivs;
  evalShLocalDerivative(coord1,coord2,coord3,shLocalDerivs);

  // Compute Jacobian Matrix
  femDoubleMat jacMat;
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,jacMat);

  // Invert Jacobian Matrix
  femDoubleMat invJacMat;
  double detJ;
  femUtils::invert3x3Matrix(jacMat,invJacMat,detJ);

  // Obtain Global SF Derivatives
  femDoubleVec temp;
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      shapeDeriv[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<kDims;loopC++){
        shapeDeriv[loopA][loopB] += invJacMat[loopB][loopC] * shLocalDerivs[loopA][loopC];
      }
    }
  }
}

// EVAL JACOBIAN DETERMINANT
double femHexa8::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  throw femException("Not Implemented.\n");
}

// EVAL JACOBIAN DETERMINANT
void femHexa8::evalJacobianMatrix(double coord1, double coord2, double coord3, femDoubleMat shDerivs){
  throw femException("Not Implemented.\n");
}

// EVAL INTEGRATION RULE
femIntegrationRule* femHexa8::evalIntegrationRule(intRuleType type){
  throw femException("Not Implemented.\n");
}

