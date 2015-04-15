#include "femNode.h"
#include "femElement.h"
#include "femFace.h"

#include "femException.h"
#include "femConstants.h"
#include "femUtils.h"

// Constructor
femElement::femElement(int number, int prop, int totalNodes, int* connections)
{
  // Assign Initial Values
  elementNumber = number;
  propertyNumber = prop;
  for(int loopA=0;loopA<totalNodes;loopA++){
    elementConnections.push_back(connections[loopA]);
  }
}

// Copy Constructor
femElement::femElement(const femElement* other){
    // Numbers
    elementNumber = other->elementNumber;
    propertyNumber = other->propertyNumber;
    // Element Connectioes
    for(unsigned int loopA=0;loopA<other->elementConnections.size();loopA++){
      elementConnections.push_back(other->elementConnections[loopA]);
    }
    // Element Faces
    for(unsigned int loopA=0;loopA<other->elementFaces.size();loopA++){
      elementFaces.push_back(other->elementFaces[loopA]);
    }
}

// ==========
// DISTRUCTOR
// ==========
femElement::~femElement(){
}

// ================================
// VIRTUAL FUNCTIONS FOR MAIN CLASS
// ================================
void femElement::fixConnectivities(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
bool femElement::isNodeInsideElement(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
void femElement::EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){
  throw femException("Not Implemented.\n");
}
bool femElement::is2D(){
  throw femException("Not Implemented.\n");
}
double femElement::EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
double femElement::EvalMixProduct(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}

// =========================================
// CHECK THE MINIMUM JACOBIAN OF THE ELEMENT
// =========================================
double femElement::checkMinDetJ(std::vector<femNode*> &nodeList, femIntegrationRule* rule){

  // Get Gauss Coords and Weightds
  femDoubleMat intCoords = rule->getCoords(numberOfNodes);

  // Loop through the Gauss Points
  double minDetJ = std::numeric_limits<double>::max();
  double detJ = 0.0;
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){
    // Eval the Jacobian Determinant at the current location
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);
    if(detJ < minDetJ){
      minDetJ = detJ;
    }
  }
  // Eval the Jacobian Determinant at the origin
  // Only for Hexaedral Elements
  if(numberOfNodes == 8){
    detJ = evalJacobian(nodeList,0.0,0.0,0.0);
    if(detJ < minDetJ){
      minDetJ = detJ;
    }
  }

  // Return Minimum Jacobian Determinant
  return minDetJ;
}


// =====================================
// VIRTUAL FUNCTIONS FOR FINITE ELEMENTS
// =====================================
void femElement::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femElement::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}

// EVAL JACOBIAN MATRIX
void evalJacobianMatrixLocal(int numberOfNodes, femDoubleMat elNodeCoords, femDoubleMat shLocalDerivs, femDoubleMat &jacMat){
  jacMat.resize(kDims);
  for(int loopA=0;loopA<kDims;loopA++){
    jacMat[loopA].resize(kDims);
  }
  // Assemble Jacobian Matrix
  for(int loopA=0;loopA<kDims;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      jacMat[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        jacMat[loopA][loopB] += shLocalDerivs[loopC][loopB] * elNodeCoords[loopC][loopA];
      }
    }
  }
}

// ====================
// EVAL JACOBIAN MATRIX
// ====================
void femElement::evalJacobianMatrix(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &jacMat){

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
  evalLocalShapeFunctionDerivative(nodeList,coord1,coord2,coord3,shLocalDerivs);

  // Compute Jacobian Matrix
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,jacMat);
}

// ===========================
// COMMON ELEMENT CONSTRUCTION
// ===========================
double femElement::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){

  // Compute Node Coordinate Vector
  int currNode = 0;
  femDoubleMat elNodeCoords;
  elNodeCoords.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
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
  evalLocalShapeFunctionDerivative(nodeList,coord1,coord2,coord3,shLocalDerivs);

  // Compute Jacobian Matrix
  femDoubleMat jacMat;
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,jacMat);

  // Print Jac Matrix
  //printf("JACOBIAN MATRIX\n");
  //for(int loopA=0;loopA<kDims;loopA++){
  //  for(int loopB=0;loopB<kDims;loopB++){
  //    printf("%e ",jacMat[loopA][loopB]);
  //  }
  //  printf("\n");
  //}
  //printf("\n");

  // Invert Jacobian Matrix
  femDoubleMat invJacMat;
  double detJ;
  femUtils::invert3x3Matrix(jacMat,invJacMat,detJ);

  // Return
  return detJ;
}

// ======================================
// EVAL GLOBAL SHAPE FUNCTION DERIVATIVES
// ======================================
void femElement::evalGlobalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &globShDeriv){

  // Compute Node Coordinate Vector
  int currNode = 0;
  femDoubleMat elNodeCoords;
  elNodeCoords.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
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
  evalLocalShapeFunctionDerivative(nodeList,coord1,coord2,coord3,shLocalDerivs);

  //printf("Local Shape Functions\n");
  //for(int loopA=0;loopA<numberOfNodes;loopA++){
  //    printf("Node %d, dx: %f, dy: %f, dz: %f\n",loopA,shLocalDerivs[loopA][0],shLocalDerivs[loopA][1],shLocalDerivs[loopA][2]);
  //}
  //printf("\n");

  // Compute Jacobian Matrix
  femDoubleMat jacMat;
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,jacMat);

  // Invert Jacobian Matrix
  femDoubleMat invJacMat;
  double detJ;
  femUtils::invert3x3Matrix(jacMat,invJacMat,detJ);

  //printf("Jacobian\n");
  //for(int loopA=0;loopA<kDims;loopA++){
  //    printf("%f %f %f\n",jacMat[loopA][0],jacMat[loopA][1],jacMat[loopA][2]);
  //}
  //printf("\n");


  // Allocate Global Shape derivatives
  globShDeriv.clear();
  globShDeriv.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    globShDeriv[loopA].resize(kDims);
  }

  // Obtain Global SF Derivatives
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      globShDeriv[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<kDims;loopC++){
        // Transpose
        globShDeriv[loopA][loopB] += invJacMat[loopB][loopC] * shLocalDerivs[loopA][loopC];
      }
    }
  }
}

// ====================
// CREATE BOUNDARY FACE
// ====================
femElement* createBoudaryFace(int faceID){

}


// =====================
// Eval Element Centroid
// =====================
void femElement::evalElementCentroid(std::vector<femNode*> &nodeList, double* centroid){
  // Initialize Centroid
  centroid[0] = 0.0;  centroid[1] = 0.0;
  centroid[2] = 0.0;
  int currNode = 0;
  int totNodes = elementConnections.size();
  for(int loopA=0;loopA<totNodes;loopA++){
    currNode = elementConnections[loopA];
    centroid[0] = centroid[0] + nodeList[currNode]->coords[0];
    centroid[1] = centroid[1] + nodeList[currNode]->coords[1];
    centroid[2] = centroid[2] + nodeList[currNode]->coords[2];
  }
  centroid[0] = centroid[0]/double(totNodes);
  centroid[1] = centroid[1]/double(totNodes);
  centroid[2] = centroid[2]/double(totNodes);
}

// ====================================================
// Eval Distance between the centroid and a given point
// ====================================================
double femElement::evalPointToElementDistance(double* pointCoords, std::vector<femNode*> &nodeList){
  double centroid[3] = {0.0};
  evalElementCentroid(nodeList,centroid);
  double dist = 0.0;
  for(int loopA=0;loopA<3;loopA++){
    dist += (centroid[loopA]-pointCoords[loopA])*(centroid[loopA]-pointCoords[loopA]);
  }
  // Return
  return sqrt(dist);
}

// ============================================================
// Eval TETRA10 Volume Coordinates: STRAIGHT SIDE APPROXIMATION
// ============================================================
void femTetra10::EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){
  int* connections= new int[4];

  // Copy the first four connctions !!! Complete...
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    connections[loopA] = elementConnections[loopA];
  }

  // Create a temporay Tet4 Element with the First 4 nodes
  femTetra4* tet4 = new femTetra4(1,1,kTetra4Nodes,connections);
  // Compute Final Volume Coordinates
  double tet4VolCoords[kTetra4Nodes] = {0.0};
  tet4->EvalVolumeCoordinates(dispFactor,pointCoords,nodeList,tet4VolCoords);

  // Compute The Area Coordinates for the Full Quadratic Tetrahedron
  volCoords[0] = tet4VolCoords[0]*(2.0*tet4VolCoords[0]-1.0);
  volCoords[1] = tet4VolCoords[1]*(2.0*tet4VolCoords[1]-1.0);
  volCoords[2] = tet4VolCoords[2]*(2.0*tet4VolCoords[2]-1.0);
  volCoords[3] = tet4VolCoords[3]*(2.0*tet4VolCoords[3]-1.0);
  volCoords[4] = 4.0*tet4VolCoords[0]*tet4VolCoords[1];
  volCoords[5] = 4.0*tet4VolCoords[1]*tet4VolCoords[2];
  volCoords[6] = 4.0*tet4VolCoords[2]*tet4VolCoords[0];
  volCoords[7] = 4.0*tet4VolCoords[0]*tet4VolCoords[3];
  volCoords[8] = 4.0*tet4VolCoords[1]*tet4VolCoords[3];
  volCoords[9] = 4.0*tet4VolCoords[2]*tet4VolCoords[3];

  // Delete Temporary Tet4 element
  delete [] connections;
  delete tet4;
}





// =================================
// Interpolate Element Displacements
// =================================
void femElement::InterpolateElementDisplacements(double dispFactor, double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps){

  double* volCoords = new double[int(elementConnections.size())];
  // Eval Shape Function Values for Current Node
  EvalVolumeCoordinates(dispFactor,nodeCoords,nodeList,volCoords);

  // Loop through the nodes
  for(int loopA=0;loopA<kNodeDofs;loopA++){
    intDisps[loopA] = 0.0;
  }

  // Loop through the nodes
  int currNode = 0;
  double* currNodeDisps = nullptr;
  for(unsigned int loopA=0;loopA<elementConnections.size();loopA++){
    // Get Node Number
    currNode = elementConnections[loopA];
    // Get Node Displacements
    currNodeDisps = nodeList[currNode]->getNodeDisplacements();
    // Loop through the node dofs
    for(int loopB=0;loopB<kNodeDofs;loopB++){
      intDisps[loopB] += currNodeDisps[loopB]*volCoords[loopA];
    }
  }
}

// =========================================
// FINITE ELEMENT ROUTINES FOR TET10 ELEMENT
// =========================================
void femTetra10::fixConnectivities(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
void femTetra10::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femTetra10::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,
                                                  femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}

// =============================================
// CAREFUL: VALID ONLY FOR STRAIGHT SIDE TETRA10
// =============================================
bool femTetra10::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  // Init Result
  bool isInside = true;

  // Copy the first four connctions
  int* connections= new int[4];
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    connections[loopA] = this->elementConnections[loopA];
  }

  // Create a temporay Tet4 Element with the First 4 nodes
  femTetra4* tet4 = new femTetra4(1,1,kTetra4Nodes,connections);
  // Compute Final Volume Coordinates
  double tet4VolCoords[kTetra4Nodes] = {0.0};
  tet4->EvalVolumeCoordinates(dispFactor,pointCoords,nodeList,tet4VolCoords);

  // Compute Result
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    isInside = (isInside)&&(tet4VolCoords[loopA] > -1.0e-2);
  }

  // Check if the volume coordinates sum up to one
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    sum += tet4VolCoords[loopA];
  }
  // NO: I use it to evaluate is a node is inside
  //if (fabs(sum-1.0)>kMathZero){
  //  throw femException("Error: Internal. In Tet10 internal Check, Shape Function don't sum up to one");
  //}

  // Delete Temporary Tet4 element
  delete [] connections;
  delete tet4;

  // Return
  return isInside;
}

// ============================================
// Check Criterion for Random Walking algorithm
// ============================================
void femElement::CheckRandomWalkingCriterion(int localFaceID, double* nodeCoords,
                                             std::vector<femFace*> &faceList,std::vector<femNode*> &nodeList,
                                             bool &isOnFace, bool &isOnOppositeSide){
  // Intialize Face Node Coords
  double faceNodeCoords[3][3];
  double normal[3] = {0.0};
  double centroid[3] = {0.0};

  // Eval element Centroid
  evalElementCentroid(nodeList,centroid);

  // Get global face number
  int faceID = elementFaces[localFaceID];

  // Get node Coordinates from face
  int currNode = 0;
  for(int loopA=0;loopA<3;loopA++){
    currNode = faceList[faceID]->faceNodes[loopA];
    faceNodeCoords[0][loopA] = nodeList[currNode]->coords[0];
    faceNodeCoords[1][loopA] = nodeList[currNode]->coords[1];
    faceNodeCoords[2][loopA] = nodeList[currNode]->coords[2];
  }

  // Get Two Vectors
  double vec1[3] = {0.0};
  double vec2[3] = {0.0};
  for(int loopA=0;loopA<3;loopA++){
    vec1[loopA] = faceNodeCoords[loopA][1] - faceNodeCoords[loopA][0];
    vec2[loopA] = faceNodeCoords[loopA][2] - faceNodeCoords[loopA][0];
  }

  // Make external Product
  femUtils::Do3DExternalProduct(vec1,vec2,normal);

  // Check the Coordinate along the normal for the point and centroid
  double pointDist = femUtils::Do3DInternalProduct(nodeCoords,normal);
  double elDist = femUtils::Do3DInternalProduct(centroid,normal);

  // First Check
  isOnFace = (fabs(pointDist)<kMathZero);

  // Second Check
  isOnOppositeSide = ((pointDist/elDist) < 0.0);
}

// =================================
// Get Adjacent Element through Face
// =================================
int femElement::getAdjacentElement(int localFaceID, std::vector<femFace*> &faceList){
  int currNumber = elementNumber;
  int faceID = elementFaces[localFaceID];
  int numElements = faceList[faceID]->faceElements.size();
  if(numElements == 1){
    return -1;
  }else{
    int el1 = faceList[faceID]->faceElements[0];
    int el2 = faceList[faceID]->faceElements[1];
    if(el1 == currNumber){
      return el2;
    }else if(el2 == currNumber){
      return el1;
    }else{
      throw femException("Internal: Invalid Element connections in face.\n");
    }
  }
}

// ==========================
// Get Bounding Box Node List
// ==========================
void femElement::CreateBoundingBoxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &boxNodeList){
  double elBox[6] = {0.0};
  double coord[3] = {0.0};
  double currDisps[6] = {0.0};
  int currNode = 0;
  // Create Element Box
  for(unsigned int loopA=0;loopA<elementConnections.size();loopA++){
    // Get current node
    currNode = elementConnections[loopA];
    // Min X
    if (nodeList[currNode]->coords[0]<elBox[0]){
      elBox[0] = nodeList[currNode]->coords[0];
    }
    // Max X
    if (nodeList[currNode]->coords[0]>elBox[1]){
      elBox[1] = nodeList[currNode]->coords[0];
    }
    // Min Y
    if (nodeList[currNode]->coords[1]<elBox[2]){
      elBox[2] = nodeList[currNode]->coords[1];
    }
    // Max Y
    if (nodeList[currNode]->coords[1]>elBox[3]){
      elBox[3] = nodeList[currNode]->coords[1];
    }
    // Min Z
    if (nodeList[currNode]->coords[2]<elBox[4]){
      elBox[4] = nodeList[currNode]->coords[2];
    }
    // Max Z
    if (nodeList[currNode]->coords[2]>elBox[5]){
      elBox[5] = nodeList[currNode]->coords[2];
    }
  }
  // Write Coords for every point
  for(int loopA=0;loopA<8;loopA++){
    switch(loopA){
      case 0:
        coord[0] = elBox[0];
        coord[1] = elBox[2];
        coord[2] = elBox[4];
        break;
      case 1:
        coord[0] = elBox[1];
        coord[1] = elBox[2];
        coord[2] = elBox[4];
        break;
      case 2:
        coord[0] = elBox[0];
        coord[1] = elBox[3];
        coord[2] = elBox[4];
        break;
      case 3:
        coord[0] = elBox[1];
        coord[1] = elBox[3];
        coord[2] = elBox[4];
        break;
      case 4:
        coord[0] = elBox[0];
        coord[1] = elBox[2];
        coord[2] = elBox[5];
        break;
      case 5:
        coord[0] = elBox[1];
        coord[1] = elBox[2];
        coord[2] = elBox[5];
        break;
      case 6:
        coord[0] = elBox[0];
        coord[1] = elBox[3];
        coord[2] = elBox[5];
        break;
      case 7:
        coord[0] = elBox[1];
        coord[1] = elBox[3];
        coord[2] = elBox[5];
        break;
    }
    // Add to Node List
    femNode* newNode = new femNode(loopA,coord,currDisps);
    boxNodeList.push_back(newNode);
  }
}

// ==============================
// Get Max-Min Bounding Box Nodes
// ==============================
void femElement::CreateMinMaxNodeList(std::vector<femNode*> &nodeList,std::vector<femNode*> &minMaxNodeList){
  double elBox[6] = {0.0};
  elBox[0] = std::numeric_limits<double>::max();
  elBox[1] = -std::numeric_limits<double>::max();
  elBox[2] = std::numeric_limits<double>::max();
  elBox[3] = -std::numeric_limits<double>::max();
  elBox[4] = std::numeric_limits<double>::max();
  elBox[5] = -std::numeric_limits<double>::max();
  double coord[3] = {0.0};
  double currDisps[6] = {0.0};
  int currNode = 0;
  // Create Element Box
  for(unsigned int loopA=0;loopA<elementConnections.size();loopA++){
    // Get current node
    currNode = elementConnections[loopA];
    // Min X
    if (nodeList[currNode]->coords[0]<elBox[0]){
      elBox[0] = nodeList[currNode]->coords[0];
    }
    // Max X
    if (nodeList[currNode]->coords[0]>elBox[1]){
      elBox[1] = nodeList[currNode]->coords[0];
    }
    // Min Y
    if (nodeList[currNode]->coords[1]<elBox[2]){
      elBox[2] = nodeList[currNode]->coords[1];
    }
    // Max Y
    if (nodeList[currNode]->coords[1]>elBox[3]){
      elBox[3] = nodeList[currNode]->coords[1];
    }
    // Min Z
    if (nodeList[currNode]->coords[2]<elBox[4]){
      elBox[4] = nodeList[currNode]->coords[2];
    }
    // Max Z
    if (nodeList[currNode]->coords[2]>elBox[5]){
      elBox[5] = nodeList[currNode]->coords[2];
    }
  }
  // Write Coords for every point
  for(int loopA=0;loopA<2;loopA++){
    switch(loopA){
      case 0:
        coord[0] = elBox[0];
        coord[1] = elBox[2];
        coord[2] = elBox[4];
        break;
      case 1:
        coord[0] = elBox[1];
        coord[1] = elBox[3];
        coord[2] = elBox[5];
        break;
    }
    // Add to Node List
    femNode* newNode = new femNode(loopA,coord,currDisps);
    minMaxNodeList.push_back(newNode);
  }
}

// =====================================
// Eval Mixed Product of Single Elements
// =====================================
double femTetra4::EvalMixProduct(std::vector<femNode*> &nodeList){
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
  vecC[0] = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  vecC[1] = nodeList[elementConnections[3]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  vecC[2] = nodeList[elementConnections[3]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  // Get external product
  femUtils::Do3DExternalProduct(vecA,vecB,vec1);
  // Get Internal Product
  currProd = femUtils::Do3DInternalProduct(vec1,vecC);
  return currProd;
}

// TO IMPLEMENT
double femTetra10::EvalMixProduct(std::vector<femNode*> &nodeList){
  throw femException("Internal: femTetra10::EvalMixProduct not Implemented");
  return 0.0;
}

// TO IMPLEMENT
double femTri3::EvalMixProduct(std::vector<femNode*> &nodeList){
  throw femException("Internal: femTria3::EvalMixProduct not Implemented");
  return 0.0;
}

// ================================
// FINITE ELEMENT ROUTINES FOR TRI3
// ================================
// =====================================
// VIRTUAL FUNCTIONS FOR FINITE ELEMENTS
// =====================================
void femTri3::fixConnectivities(std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
void femTri3::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  shapeFunction.clear();
  shapeFunction.push_back(coord1);
  shapeFunction.push_back(coord2);
  shapeFunction.push_back(1.0-coord1-coord2);
}
void femTri3::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  femDoubleVec temp;
  // Third Component Set to Zero: CHECK!!!
  temp.push_back(1.0);
  temp.push_back(0.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
  temp.push_back(0.0);
  temp.push_back(1.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
  temp.push_back(-1.0);
  temp.push_back(-1.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
}

// ====================================================
// ASSEMBLE MASS MATRIX FOR ADVECTION DIFFUSION PROBLEM
// ====================================================
void femElement::assembleMass(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double** massMat){
  femDoubleVec shapeFunction;
  femDoubleMat shapeFunctionDerivs;
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  int globNode = 0;
  double currProd = 0.0;

  // Get Gauss Points and Weights
  intCoords = rule.getCoords(numberOfNodes);
  intWeights = rule.getWeights(numberOfNodes);

  for(int loopA=0;loopA<rule.getTotGP(numberOfNodes);loopA++){

    // Eval Shape Functions
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);
    // Eval Shape Function Derivatives
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunctionDerivs);

    // Loop Through the Nodes
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        // Eval Product with Velocity
        currProd = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          globNode = elementConnections[loopB];
          currProd += nodeVelocities[globNode][loopD] * shapeFunctionDerivs[loopB][loopD];
        }
        massMat[loopB][loopC] += shapeFunction[loopB]*shapeFunction[loopC] + tauSUPG[loopA] * currProd * shapeFunction[loopC] * intWeights[loopA];
      }
    }
  }
}

// =========================================================
// ASSEMBLE STIFFNESS MATRIX FOR ADVECTION DIFFUSION PROBLEM
// =========================================================
void femElement::assembleStiffness(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double diffusivity, femDoubleMat &stiffnessMat){
  femDoubleVec shapeFunction;
  femDoubleMat shapeFunctionDerivs;
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  int globNode = 0;
  double uNablaNB = 0.0;
  double uNablaNA = 0.0;
  double NablaNANablaNB = 0.0;
  double firstTerm = 0.0;
  double secondTerm = 0.0;
  double tauTerm = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule.getCoords(numberOfNodes);
  intWeights = rule.getWeights(numberOfNodes);

  for(int loopA=0;loopA<rule.getTotGP(numberOfNodes);loopA++){

    // Eval Shape Functions
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);
    // Eval Shape Function Derivatives
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunctionDerivs);
    // Loop Through the Nodes
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        // Eval Product with Velocity
        uNablaNB = 0.0;
        uNablaNA = 0.0;
        NablaNANablaNB = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          globNode = elementConnections[loopB];
          uNablaNA += nodeVelocities[globNode][loopD] * shapeFunctionDerivs[loopB][loopD];
          uNablaNB += nodeVelocities[globNode][loopD] * shapeFunctionDerivs[loopC][loopD];
          NablaNANablaNB += shapeFunctionDerivs[loopB][loopD] * shapeFunctionDerivs[loopC][loopD];
        }
        firstTerm = shapeFunction[loopA] * uNablaNB;
        secondTerm = diffusivity * NablaNANablaNB;
        tauTerm = tauSUPG[loopA] * uNablaNA * uNablaNB;
        // Weights???
        stiffnessMat[loopB][loopC] += firstTerm + secondTerm + tauTerm;
      }
    }
  }
}

// =================================
// POISSON PROBLEM - ASSEMBLE MATRIX
// =================================
void femElement::formPoissonMatrix(std::vector<femNode*> nodeList,femIntegrationRule* rule,femDoubleVec diffusivity,femDoubleMat &elMat){

  // CLEAR AND ALLOCATE MATRIX
  elMat.clear();
  elMat.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elMat[loopA].resize(numberOfNodes);
  }
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      elMat[loopA][loopB] = 0.0;
    }
  }

  //printf("Diff: %e %e %e\n",diffusivity[0],diffusivity[1],diffusivity[2]);

  // INIT SHAPE DERIVATIVE MATRIX
  femDoubleMat shapeDeriv;

  // GAUSS POINTS LOOP
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;
  double currStiff = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes);
  intWeights = rule->getWeights(numberOfNodes);

  //printf("Gauss Points: %d\n",intCoords.size());
  //for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){
  //  printf("GP: %d, gpx: %f, gpy: %f, gpz: %f\n",loopA,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);
  //}
  //printf("\n");


  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){

    // Eval Current Shape Derivatives Matrix
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeDeriv);

    // Print Global Derivatives
    //printf("Global Derivatives\n");
    //for(int loopA=0;loopA<numberOfNodes;loopA++){
    //  printf("Node %d, dx: %f, dy: %f, dz: %f\n",loopA,shapeDeriv[loopA][0],shapeDeriv[loopA][1],shapeDeriv[loopA][2]);
    //}
    //printf("\n");

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    // Eval Resulting Matrix
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        currStiff = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          currStiff += diffusivity[loopD] * shapeDeriv[loopB][loopD] * shapeDeriv[loopC][loopD];
        }
        // Assemble Gauss Point Contribution
        elMat[loopB][loopC] += currStiff * detJ * intWeights[loopA];
      }
    }
  }
}

// =================================
// POISSON PROBLEM - ASSEMBLE SOURCE
// =================================
void femElement::formPoissonSource(std::vector<femNode*> nodeList,femIntegrationRule* rule,double sourceValue,femDoubleVec &elSourceVec){

  // Shape Function Values
  femDoubleVec shapeFunction;

  // Init Source Array
  elSourceVec.clear();
  elSourceVec.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elSourceVec[loopA] = 0.0;
  }

  // Gauss Point Loop
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes);
  intWeights = rule->getWeights(numberOfNodes);

  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){

    // Eval Shape Function At the current GP
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    // Eval Resulting Matrix
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      elSourceVec[loopB] += shapeFunction[loopB] * sourceValue * detJ * intWeights[loopA];
    }
  }
}

// =================================
// POISSON PROBLEM - ASSEMBLE SOURCE
// =================================
void femElement::formPoissonNeumannBC(){
}

// ======================================
// INTEGRATE NODAL QUANTITY OVER FEM MESH
// ======================================
double femElement::integrateNodalVector(std::vector<femNode*> nodeList,femIntegrationRule* rule,femDoubleVec nodeVec){
  // Get Gauss Coords and Weightds
  femDoubleMat intCoords = rule->getCoords(numberOfNodes);
  femDoubleVec intWeights = rule->getWeights(numberOfNodes);

  // Eval Shape Functions
  femDoubleVec shapeFunction;

  // Loop on the Gauss Point locations
  double result = 0.0;
  double detJ = 0.0;
  double currGPValue = 0.0;
  int currNode = 0;
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){
    // Eval the shape Functions at the current Integration Point
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);
    // Eval the Jacobian Determinant at the current location
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);
    // Eval Current Function Value
    currGPValue = 0.0;
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      currNode = elementConnections[loopB];
      currGPValue += shapeFunction[loopB] * nodeVec[currNode];
    }
    // Sum Gauss Point Contribution
    result += currGPValue * intWeights[loopA] * detJ;
  }
  //printf("GP %d, Value %f, Weight %f, Result: %f\n",rule->getTotGP(numberOfNodes),currGPValue,intWeights[0],result);
  return result;
}

// ASSEMBLE ADVECTION DIFFUSION EQUATION
void femElement::formAdvDiffLHS(std::vector<femNode*> nodeList,femIntegrationRule* rule,femDoubleVec diffusivity,femDoubleVec velocity,femDoubleMat &elMat){

  // GET SCALAR Diffusivity and Velocity
  scalarDiff = diffusivity[0];
  scalarVel = velocity[0];

  // CLEAR AND ALLOCATE MATRIX
  elMat.clear();
  elMat.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elMat[loopA].resize(numberOfNodes);
  }
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      elMat[loopA][loopB] = 0.0;
    }
  }

  // INIT SHAPE DERIVATIVE MATRIX
  femDoubleMat shapeDeriv;

  // GAUSS POINTS LOOP
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;
  double currStiff = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes);
  intWeights = rule->getWeights(numberOfNodes);

  // Gauss Point Loop
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){

    // Eval Shape Function
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);

    // Eval Current Shape Derivatives Matrix
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeDeriv);

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    // Eval Resulting Matrix
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        currStiff = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          currStiff +=
          // Diffusion
          scalarDiff * shapeDeriv[loopB][loopD] * shapeDeriv[loopC][loopD] +
          // Advection
          scalarVel * shapeFunction[loopB] * shapeDeriv[loopC][loopD];
        }
        // Assemble Gauss Point Contribution
        elMat[loopB][loopC] += currStiff * detJ * intWeights[loopA];
      }
    }
  }
}

void femElement::formAdvDiffRHS(std::vector<femNode*> nodeList,femIntegrationRule* rule,double sourceValue,femDoubleVec &elSourceVec){

  // Shape Function Values
  femDoubleVec shapeFunction;

  // Init Source Array
  elSourceVec.clear();
  elSourceVec.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elSourceVec[loopA] = 0.0;
  }

  // Gauss Point Loop
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes);
  intWeights = rule->getWeights(numberOfNodes);

  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){

    // Eval Shape Function At the current GP
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    // Eval Resulting Matrix
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      elSourceVec[loopB] += shapeFunction[loopB] * sourceValue * detJ * intWeights[loopA];
    }
  }
}

