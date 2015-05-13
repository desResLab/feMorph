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
  femDoubleMat intCoords = rule->getCoords(numberOfNodes,dims);

  // Loop through the Gauss Points
  double minDetJ = std::numeric_limits<double>::max();
  double detJ = 0.0;
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){
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

// ====================
// EVAL JACOBIAN MATRIX
// ====================
void evalJacobianMatrixLocal(int numberOfNodes, femDoubleMat elNodeCoords, femDoubleMat shLocalDerivs, elDim dims, femDoubleMat &jacMat){
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
  // Put one on the diagonal for 1D and 2D problems
  if(dims == d1){
    jacMat[1][1] = 1.0;
    jacMat[2][2] = 1.0;
  }else if(dims == d2){
    jacMat[2][2] = 1.0;
  }
}

// ====================
// EVAL JACOBIAN MATRIX
// ====================
void femElement::evalJacobianMatrix(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, elDim dims, femDoubleMat &jacMat){

  // Compute Node Coordinate Vector
  int currNode = 0;
  femDoubleMat elNodeCoords;
  elNodeCoords.resize(8);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elNodeCoords[loopA].resize(kDims);
  }
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    currNode = elementConnections[loopA];
    for(int loopB=0;loopB<kDims;loopB++){
      elNodeCoords[loopA][loopB] = nodeList[currNode]->coords[loopB];
      printf("Curr Coord: %f",nodeList[currNode]->coords[loopB]);
    }
  }

  printf("CIAO\n");

  // Compute Local Derivatives
  femDoubleMat shLocalDerivs;
  evalLocalShapeFunctionDerivative(nodeList,coord1,coord2,coord3,shLocalDerivs);

  // Compute Jacobian Matrix
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,dims,jacMat);
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
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,dims,jacMat);

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

// =============================
// EVAL GEOMETRIC ELEMENT MATRIX
// =============================
void femElement::evalGeometricMatrix(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,femDoubleMat& elGeomMat){
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
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,dims,jacMat);

  // Invert Jacobian Matrix
  femDoubleMat invJacMat;
  double detJ;
  femUtils::invert3x3Matrix(jacMat,invJacMat,detJ);

  // Allocate Geometric Matrix
  elGeomMat.resize(kDims);
  for(int loopA=0;loopA<kDims;loopA++){
    elGeomMat[loopA].resize(kDims);
  }

  // Compute Geometric Matrix
  for(int loopA=0;loopA<kDims;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      elGeomMat[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<kDims;loopC++){
        // CAREFULL !!!
        elGeomMat[loopA][loopB] += invJacMat[loopC][loopA] * invJacMat[loopC][loopB];
      }
    }
  }
}

// ==========================
// EVAL QUADRATIC FORM OF TAU
// ==========================
double evalQuadraticTau(femDoubleMat elGeomMat,femDoubleVec velocity,femDoubleVec diffusivity){
  // Select Constant
  double cConst = 1.0;
  double isoDiff = sqrt(diffusivity[0]*diffusivity[0] + diffusivity[1]*diffusivity[1] + diffusivity[2]*diffusivity[2]);

  // Contribution from Advection
  double firstTau = 0.0;
  for(int loopA=0;loopA<kDims;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      firstTau += velocity[loopA] * elGeomMat[loopA][loopB] * velocity[loopB];
    }
  }

  // Contribution from Diffusion
  double secondTau = 0.0;
  for(int loopA=0;loopA<kDims;loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      secondTau += elGeomMat[loopA][loopB] * elGeomMat[loopA][loopB];
    }
  }
  secondTau *= cConst * cConst * isoDiff * isoDiff;

  //printf("Tau 1 %f, Tau 2 %f\n",firstTau,secondTau);

  // Return the norm
  return 1.0/sqrt(firstTau + secondTau);
}

// ======================================
// EVAL GLOBAL SHAPE FUNCTION DERIVATIVES
// ======================================
void femElement::evalGlobalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &globShDeriv){

  // Flag to print shape functions
  bool printSF = false;
  bool printJAC = false;

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

  // Print Shape Functions if requested
  if(printSF){
    printf("Local Shape Functions\n");
    for(int loopA=0;loopA<numberOfNodes;loopA++){
        printf("Node %d, dx: %f, dy: %f, dz: %f\n",loopA,shLocalDerivs[loopA][0],shLocalDerivs[loopA][1],shLocalDerivs[loopA][2]);
    }
    printf("\n");
  }

  // Compute Jacobian Matrix
  femDoubleMat jacMat;
  evalJacobianMatrixLocal(numberOfNodes,elNodeCoords,shLocalDerivs,dims,jacMat);

  // Invert Jacobian Matrix
  femDoubleMat invJacMat;
  double detJ;
  femUtils::invert3x3Matrix(jacMat,invJacMat,detJ);

  if(printJAC){
    printf("Jacobian\n");
    for(int loopA=0;loopA<kDims;loopA++){
      printf("%12.3f %12.3f %12.3f\n",jacMat[loopA][0],jacMat[loopA][1],jacMat[loopA][2]);
    }
    printf("\n");
  }

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
  intCoords = rule.getCoords(numberOfNodes,dims);
  intWeights = rule.getWeights(numberOfNodes,dims);

  for(int loopA=0;loopA<rule.getTotGP(numberOfNodes,dims);loopA++){

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
  intCoords = rule.getCoords(numberOfNodes,dims);
  intWeights = rule.getWeights(numberOfNodes,dims);

  for(int loopA=0;loopA<rule.getTotGP(numberOfNodes,dims);loopA++){

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
  intCoords = rule->getCoords(numberOfNodes,dims);
  intWeights = rule->getWeights(numberOfNodes,dims);

  //printf("Gauss Points: %d\n",intCoords.size());
  //for(int loopA=0;loopA<rule->getTotGP(numberOfNodes);loopA++){
  //  printf("GP: %d, gpx: %f, gpy: %f, gpz: %f\n",loopA,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);
  //}
  //printf("\n");


  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){

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
  intCoords = rule->getCoords(numberOfNodes,dims);
  intWeights = rule->getWeights(numberOfNodes,dims);

  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){

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
  femDoubleMat intCoords = rule->getCoords(numberOfNodes,dims);
  femDoubleVec intWeights = rule->getWeights(numberOfNodes,dims);

  // Eval Shape Functions
  femDoubleVec shapeFunction;

  // Loop on the Gauss Point locations
  double result = 0.0;
  double detJ = 0.0;
  double currGPValue = 0.0;
  int currNode = 0;
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){
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

// ================================
// ASSEMBLE ADVECTION DIFFUSION LHS
// ================================
void femElement::formAdvDiffLHS(std::vector<femNode*> nodeList,
                                femIntegrationRule* rule,
                                femDoubleVec diffusivity,femDoubleVec velocity,
                                int schemeType,
                                femDoubleMat &elMat){

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
  femDoubleVec shapeFunction;
  femDoubleMat elGeomMat;

  // GAUSS POINTS LOOP
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;
  double currTau = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes,dims);
  intWeights = rule->getWeights(numberOfNodes,dims);

  // printf("Total Gauss Points: %d\n",rule->getTotGP(numberOfNodes));

  // Gauss Point Loop
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){

    // Eval Shape Function
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);

    //printf(" Gauss Point %f %f %f\n",intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);
    //printf(" Shape Function %f %f\n",shapeFunction[0],shapeFunction[1]);

    // Eval Current Shape Derivatives Matrix
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeDeriv);

    // Eval Geometric Element Matrix
    evalGeometricMatrix(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],elGeomMat);

    // Eval Quadratic Tau
    currTau = evalQuadraticTau(elGeomMat,velocity,diffusivity);

    //printf("dN1/dx %f, dN1/dy %f, dN1/dz %f\n",shapeDeriv[0][0],shapeDeriv[0][1],shapeDeriv[0][2]);
    //printf("dN2/dx %f, dN2/dy %f, dN2/dz %f\n",shapeDeriv[1][0],shapeDeriv[1][1],shapeDeriv[1][2]);

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    //printf("detJ %f\n",detJ);

    // Eval Resulting Matrix
    double currStiff = 0.0;
    double diffTerm = 0.0;
    double advTerm = 0.0;
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        currStiff = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){          
          // SUPG
          // SUPG Diffusion: CHECK !!!
          diffTerm = diffusivity[loopD] * shapeDeriv[loopB][loopD] * shapeDeriv[loopC][loopD];
          // SUPG Advection: CHECK !!!
          advTerm = velocity[loopD] * (shapeFunction[loopB] + currTau*shapeDeriv[loopB][loopD]) * shapeDeriv[loopC][loopD];
          // Combine the two
          currStiff += diffTerm + advTerm;
        }
        // Assemble Gauss Point Contribution
        elMat[loopB][loopC] += currStiff * detJ * intWeights[loopA];
      }
    }
  }
}

// EVAL FORCING FOR ADV-DIFF EXERCISE
double evalForcing(double coord,double advVelocity,double totLength){
  double result = 0.0;
  if((coord>=0.0)&&(coord<=3.75)){
    result = (16.0*advVelocity/totLength)*(1.0 - 4.0*coord/totLength);
  }else if((coord>3.75)&&(coord<=5.0)){
    result = (16.0*advVelocity/totLength)*(-2.0 + 4.0*coord/totLength);
  }else{
    result = 0.0;
  }
  return result;
}

// ================================
// ASSEMBLE ADVECTION DIFFUSION RHS
// ================================
void femElement::formAdvDiffRHS(std::vector<femNode*> nodeList,
                                femIntegrationRule* rule,
                                double sourceValue,
                                femDoubleVec diffusivity,femDoubleVec velocity,
                                int schemeType,
                                femDoubleVec &elSourceVec){

  // CHOOSE ARTIFICIAL DIFFUSION MODEL
  //scheme: // 0-SUPG, 1-GALERKIN, 2-UD, 3-GALERKIN + EAD

  // GET SCALAR Diffusivity and Velocity
  double scalarDiff = diffusivity[0];
  double scalarVel = velocity[0];

  // Eval Element Length
  double elSize = EvalVolume(1.0,nodeList);

  // Eval Factor
  double alpha = (scalarVel*elSize)/(2.0*scalarDiff);
  double tauFactor = femUtils::coth(alpha) - 1.0/alpha;

  // Shape Function Values
  femDoubleVec shapeFunction;
  femDoubleVec shapeFunction1;
  femDoubleVec shapeFunction2;
  femDoubleMat shapeDeriv1;
  femDoubleMat shapeDeriv2;

  // Init Source Array
  elSourceVec.clear();
  elSourceVec.resize(numberOfNodes);
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    elSourceVec[loopA] = 0.0;
  }

  // Determine Type of Forcing
  bool exactForcing = false;

  if(exactForcing){
    // Use the Trapeziodal rule
    // Get end Nodes
    double node1x = nodeList[elementConnections[0]]->coords[0];
    double node2x = nodeList[elementConnections[1]]->coords[0];

    // Choose the total number of intervals
    int totalTrapzPoints = 1024;

    // Loop on trapeziodal points
    elSourceVec[0] = 0.0;
    elSourceVec[1] = 0.0;
    double intCoordX1 = 0.0;
    double intCoordX2 = 0.0;
    double normCoord1 = 0.0;
    double normCoord2 = 0.0;
    double currForcing1 = 0.0;
    double currForcing2 = 0.0;
    for(int loopA=0;loopA<totalTrapzPoints;loopA++){
      // Eval Physical Coordinates
      intCoordX1 = node1x + loopA*(node2x - node1x)/(double)totalTrapzPoints;
      intCoordX2 = node1x + (loopA + 1)*(node2x - node1x)/(double)totalTrapzPoints;
      // Eval Normalized Coords
      normCoord1 = (2.0*(intCoordX1 - node1x)/(node2x - node1x)) - 1.0;
      normCoord2 = (2.0*(intCoordX2 - node1x)/(node2x - node1x)) - 1.0;

      // Eval Shape Function At the current Trapz Point
      evalShapeFunction(nodeList,normCoord1,0.0,0.0,shapeFunction1);
      evalShapeFunction(nodeList,normCoord2,0.0,0.0,shapeFunction2);

      // Eval Current Shape Derivatives
      evalGlobalShapeFunctionDerivative(nodeList,normCoord1,0.0,0.0,shapeDeriv1);
      evalGlobalShapeFunctionDerivative(nodeList,normCoord2,0.0,0.0,shapeDeriv2);

      if(schemeType == 0){
        shapeFunction1[0] += (0.5*elSize*tauFactor)*shapeDeriv1[0][0];
        shapeFunction1[1] += (0.5*elSize*tauFactor)*shapeDeriv1[1][0];
        shapeFunction2[0] += (0.5*elSize*tauFactor)*shapeDeriv2[0][0];
        shapeFunction2[1] += (0.5*elSize*tauFactor)*shapeDeriv2[1][0];
        //printf("SUPG SH1: %f\n",(0.5*elSize*tauFactor)*shapeDeriv1[0][0]);
        //printf("SUPG SH2: %f\n",(0.5*elSize*tauFactor)*shapeDeriv1[1][0]);
        //printf("SUPG SH3: %f\n",(0.5*elSize*tauFactor)*shapeDeriv2[0][0]);
        //printf("SUPG SH4: %f\n",(0.5*elSize*tauFactor)*shapeDeriv2[1][0]);
      }

      // Eval Forcing
      currForcing1 = evalForcing(intCoordX1,scalarVel,10.0);
      currForcing2 = evalForcing(intCoordX2,scalarVel,10.0);

      // Eval RHS for element
      elSourceVec[0] += 0.5*(intCoordX2-intCoordX1)*(currForcing1*shapeFunction1[0] + currForcing2*shapeFunction2[0]);
      elSourceVec[1] += 0.5*(intCoordX2-intCoordX1)*(currForcing1*shapeFunction1[1] + currForcing2*shapeFunction2[1]);
    }

    //printf("source 0: %f\n",elSourceVec[0]);
    //printf("source 1: %f\n",elSourceVec[1]);

  }else{
    // Gauss Point Loop
    femDoubleMat intCoords;
    femDoubleVec intWeights;
    double detJ = 0.0;

    // Get Integration Coords and Weights
    intCoords = rule->getCoords(numberOfNodes,dims);
    intWeights = rule->getWeights(numberOfNodes,dims);

    for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){

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
}

// ASSEMBLE WEAK BCS FOR BOUNDARY ELEMENTS
void femElement::formWeakBC(std::vector<femNode*> nodeList,femIntegrationRule* rule,
                            femDoubleVec diffusivity,femDoubleVec velocity,femDoubleVec elNormal, double elBCValue,
                            femDoubleMat &elMat,femDoubleVec &elVec){

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




  // CHECK IF INLET OR OUTLET
  // EVAL NORM OF DIFFUSION TENSOR
  double normalVel = 0.0;
  double diffNorm = 0.0;
  for(int loopA=0;loopA<kDims;loopA++){
    normalVel += velocity[loopA] * elNormal[loopA];
    diffNorm += diffusivity[loopA];
  }
  diffNorm = sqrt(diffNorm);
  bool isInlet = (normalVel < 0.0);
  // Define Constants
  double gamma = 1.0;
  double Cb = 1.0;
  double elSize = EvalVolume(0.0,nodeList);

  // INIT SHAPE DERIVATIVE MATRIX
  femDoubleMat shapeDeriv;
  femDoubleVec shapeFunction;
  femDoubleMat elGeomMat;

  // GAUSS POINTS LOOP
  femDoubleMat intCoords;
  femDoubleVec intWeights;
  double detJ = 0.0;
  double currTau = 0.0;

  // Get Integration Coords and Weights
  intCoords = rule->getCoords(numberOfNodes,dims);
  intWeights = rule->getWeights(numberOfNodes,dims);

  // Gauss Point Loop
  for(int loopA=0;loopA<rule->getTotGP(numberOfNodes,dims);loopA++){

    // Eval Shape Function
    evalShapeFunction(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeFunction);

    // Eval Current Shape Derivatives Matrix
    evalGlobalShapeFunctionDerivative(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2],shapeDeriv);

    // Eval Determinant of the Jacobian Matrix
    detJ = evalJacobian(nodeList,intCoords[loopA][0],intCoords[loopA][1],intCoords[loopA][2]);

    // EVAL LHS CONTRIBUTION
    double currStiff = 0.0;
    double consistencyTerm = 0.0;
    double inletTerm = 0.0;
    double outletTerm = 0.0;
    double penaltyTerm = 0.0;
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        currStiff = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          // ADD WEAK BOUNDARY CONDITIONS TERMS
          // CONSISTENCY TERM
          consistencyTerm = - shapeFunction[loopB] * diffusivity[loopD] * shapeDeriv[loopC][loopD] * elNormal[loopD] +
                              shapeFunction[loopB] * velocity[loopD] * elNormal[loopD] * shapeFunction[loopB];
          if(isInlet){
            // INLET TERM
            inletTerm = - gamma * diffusivity[loopD] * shapeDeriv[loopB][loopD] * elNormal[loopD] * shapeFunction[loopC] -
                          velocity[loopD] * elNormal[loopD]  * shapeFunction[loopB] * shapeFunction[loopC];
            // OUTLET TERM
            outletTerm = 0.0;
          }else{
              // INLET TERM
              inletTerm = 0.0;
              // OUTLET TERM
              outletTerm = - gamma * diffusivity[loopD] * shapeDeriv[loopB][loopD] * elNormal[loopD] * shapeFunction[loopC];
          }
          // PENALTY TERM
          penaltyTerm = (Cb * diffNorm / elSize) * shapeFunction[loopB] * shapeFunction[loopC];

          // SUM CONTRIBUTIONS FROM ALL TERMS
          currStiff = consistencyTerm + inletTerm + outletTerm + penaltyTerm;
        }
        // Assemble Gauss Point Contribution
        elMat[loopB][loopC] += currStiff * detJ * intWeights[loopA];
      }
    }

    // EVAL RHS CONTRIBUTION
    double currRHS = 0.0;
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      currStiff = 0.0;
      for(int loopD=0;loopD<kDims;loopD++){
        // ADD WEAK BOUNDARY CONDITIONS TERMS
        if(isInlet){
          // INLET TERM
          inletTerm = - gamma * diffusivity[loopD] * shapeDeriv[loopB][loopD] * elNormal[loopD] * elBCValue -
                        velocity[loopD] * elNormal[loopD]  * shapeFunction[loopB] * elBCValue;
          // OUTLET TERM
          outletTerm = 0.0;
        }else{
          // INLET TERM
          inletTerm = 0.0;
          // OUTLET TERM
          outletTerm = - gamma * diffusivity[loopD] * shapeDeriv[loopB][loopD] * elNormal[loopD] * elBCValue;
        }
        // PENALTY TERM
        penaltyTerm = (Cb * diffNorm / elSize) * shapeFunction[loopB] * elBCValue;
        // SUM CONTRIBUTIONS FROM ALL TERMS
        currRHS = inletTerm + outletTerm + penaltyTerm;
      }
      // Assemble Gauss Point Contribution
      elVec[loopB] += currRHS * detJ * intWeights[loopA];
    }
  }
}

