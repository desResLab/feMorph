# include "femElement.h"
# include "femException.h"

void femQuad4::fixConnectivities(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
void femQuad4::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femQuad4::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}

