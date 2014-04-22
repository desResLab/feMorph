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
// Destructor
// ==========
femElement::~femElement(){
}

// ================================
// VIRTUAL FUNCTIONS FOR MAIN CLASS
// ================================
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
double femElement::EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){
  throw femException("Not Implemented.\n");
}
void femElement::assembleMass(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double** massMat){
  throw femException("Not Implemented.\n");
}
void femElement::assembleStiffness(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double diffusivity, femDoubleMat &stiffnessMat){
  throw femException("Not Implemented.\n");
}

// =====================================
// VIRTUAL FUNCTIONS FOR FINITE ELEMENTS
// =====================================
void femElement::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femElement::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}
void femElement::evalJacobianMatrix(double coord1, double coord2, double coord3, femDoubleMat shDerivs){
  throw femException("Not Implemented.\n");
}
double femElement::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  throw femException("Not Implemented.\n");
}
femIntegrationRule* femElement::evalIntegrationRule(intRuleType type){
  throw femException("Not Implemented.\n");
}

// =====================
// Eval Element Centroid
// =====================
void femElement::evalElementCentroid(std::vector<femNode*> &nodeList, double* centroid){
  // Initialize Centroid
  centroid[0] = 0.0;
  centroid[1] = 0.0;
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

// =======================================
// Assemble Tetra4 Coordinates in a Matrix
// =======================================
void femTetra4::AssembleTetCoordsMat(double dispFactor, std::vector<femNode*> &nodeList, double** coordMat){
  // Fill Matrix
  int currNode = 0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    // Get Current Node
    currNode = elementConnections[loopA];
    // Store In Matrix
    for(int loopB=0;loopB<3;loopB++){
      coordMat[loopB][loopA] = nodeList[currNode]->coords[loopB] + dispFactor * nodeList[currNode]->displacements[loopB];
    }
  }
}

// ================================
// ASSEMBLE STABILIZATION PARAMETER
// ================================
/*virtual void femTetra4::assembleStab(femOptions options, femIntegrationRule* rule,femDoubleMat &nodeVelocities, std::vector<double> &tauSUPG){
  // Loop Through the Gauss Points
  for(int loopGauss=0;loopGauss<rule->totalGP;loopGauss++){
    evalShapeFunctionDerivative(nodeList,rule->coords[0],rule->coords[1],rule->coords[2],shapeDeriv);
    // Assemble First Term
    firstTerm[loopGauss] = 0.0;
    for(int loopA=0;loopA<numberOfNodes;loopA++){
      dotProduct = 0.0;
      currNode = elementConnections[loopA];
      for(int loopB=0;loopB<kDims;loopB++){
        dotProduct += nodeVelocities[currNode][loopB] * shapeDeriv[loopA][loopB];
      }
      firstTerm[loopGauss] += fabs(dotProduct)
    }
    firstTerm[loopGauss] = (1.0/(firstTerm[loopGauss]));
    // Assemble Second Term
    secondTerm[loopGauss] = (options.deltaT/2.0);

  }
}*/

// ====================================
// Eval Tetra4 Volume from Coord Matrix
// ====================================
double EvalTetVolumeFromCoordMat(double** coordMat){
  // Allocate
  double firstVec[3] = {0.0};
  double secondVec[3] = {0.0};
  double thirdVec[3] = {0.0};

  // Get The Vectors
  // First
  firstVec[0] = coordMat[0][1] - coordMat[0][0];
  firstVec[1] = coordMat[1][1] - coordMat[1][0];
  firstVec[2] = coordMat[2][1] - coordMat[2][0];
  // Second
  secondVec[0] = coordMat[0][2] - coordMat[0][0];
  secondVec[1] = coordMat[1][2] - coordMat[1][0];
  secondVec[2] = coordMat[2][2] - coordMat[2][0];
  // Third
  thirdVec[0] = coordMat[0][3] - coordMat[0][0];
  thirdVec[1] = coordMat[1][3] - coordMat[1][0];
  thirdVec[2] = coordMat[2][3] - coordMat[2][0];

  // Eval External Product
  double auxVec[3] = {0.0};
  femUtils::Do3DExternalProduct(secondVec,thirdVec,auxVec);

  // Eval Internal Product
  return (femUtils::Do3DInternalProduct(firstVec,auxVec)/6.0);
}

// ========================
// Get External Tet Volumes
// ========================
void EvalExternalTetVolumes(double* pointCoords, double** coordMat, double* extTetVol){
  // Allocate current coordinate matrix
  double** currCoordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    currCoordMat[loopA] = new double[kTetra4Nodes];
  }
  // Eval Volume Coords
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    for(int loopB=0;loopB<kTetra4Nodes;loopB++){
      if (loopB != loopA){
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = coordMat[loopC][loopB];
        }
      }else{
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = pointCoords[loopC];
        }
      }
    }
    // Eval Volume Using CurrCoordMat
    extTetVol[loopA] = EvalTetVolumeFromCoordMat(currCoordMat);
  }
  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] currCoordMat[loopA];
  }
  delete [] currCoordMat;
}

// ===================
// Eval Element Volume
// ===================
double femTetra4::EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){
  // Put The Coordinates in a Matrix
    // Allocate current coordinate matrix
  double** coordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    coordMat[loopA] = new double[kTetra4Nodes];
  }
  AssembleTetCoordsMat(dispFactor,nodeList,coordMat);

  // Eval Tethrahedral Volume
  return EvalTetVolumeFromCoordMat(coordMat);
}



// ========================================
// Eval Volume Coordinates of Point in Tet4
// ========================================
void femTetra4::EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){

  // Put The Coordinates in a Matrix
    // Allocate current coordinate matrix
  double** coordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    coordMat[loopA] = new double[kTetra4Nodes];
  }
  AssembleTetCoordsMat(dispFactor,nodeList,coordMat);

  // Eval Tethrahedral Volume
  double TetVolume = EvalTetVolumeFromCoordMat(coordMat);

  // Eval Other Volumes
  double currVolumes[kTetra4Nodes] = {0.0};
  EvalExternalTetVolumes(pointCoords,coordMat,currVolumes);

  // Eval Shape Function Summation
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++)
  {
    sum += currVolumes[loopA];
  }
  //if (fabs(sum-1.0)>kMathZero){
  //  throw femException("Internal: Tet4 Volume Coordinates do not sum up to one.");
  //}

  // Compute Final Volume Coordinates
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    volCoords[loopA] = (currVolumes[loopA]/TetVolume);
  }

  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] coordMat[loopA];
  }
  delete [] coordMat;
}


// Check if Node is Inside
bool femTetra4::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  // Init Result
  bool isInside = true;

  // Eval Volume Coordinates
  double volCoords[kTetra4Nodes] = {0.0};
  EvalVolumeCoordinates(dispFactor,pointCoords,nodeList,volCoords);

  // Compute Result
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    isInside = (isInside)&&(volCoords[loopA] >= 0.0);
  }

  // Check if the volume coordinates sum up to one
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    sum += volCoords[loopA];
  }
  // NO: I use it to evaluate if a node is inside
  //if (fabs(sum-1.0)>kMathZero){
  //  throw femException("Error: Internal. Tet4 Shape Function don't sum up to one");
  //}
  // Return
  return isInside;
}

// ==========================================
// FINITE ELEMENT ROUTINES FOR TETRA4 ELEMENT
// ==========================================
void femTetra4::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  shapeFunction.clear();
  shapeFunction.push_back(coord1);
  shapeFunction.push_back(coord2);
  shapeFunction.push_back(coord3);
  shapeFunction.push_back(1.0-coord1-coord2-coord3);
}


void femTetra4::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,
                                            femDoubleMat &shapeDeriv){

  double x12 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x13 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[2]]->coords[0];
  double x14 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x21 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x24 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x31 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x32 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x34 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x42 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x43 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[2]]->coords[0];

  double y12 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y13 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[2]]->coords[1];
  double y14 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y21 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y24 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y31 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y32 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y34 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y42 = nodeList[elementConnections[3]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y43 = nodeList[elementConnections[3]]->coords[1] - nodeList[elementConnections[2]]->coords[1];

  double z12 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z13 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[2]]->coords[2];
  double z14 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z21 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  double z24 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z31 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  double z32 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z34 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z42 = nodeList[elementConnections[3]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z43 = nodeList[elementConnections[3]]->coords[2] - nodeList[elementConnections[2]]->coords[2];

  shapeDeriv[0][0] = y42*z32-y32*z42;
  shapeDeriv[1][0] = y31*z43-y34*z13;
  shapeDeriv[2][0] = y24*z14-y14*z24;
  shapeDeriv[3][0] = y13*z21-y12*z31;

  shapeDeriv[0][1] = x32*z42-x42*z32;
  shapeDeriv[1][1] = x43*z31-x13*z34;
  shapeDeriv[2][1] = x14*z24-x24*z14;
  shapeDeriv[3][1] = x21*z13-x31*z12;

  shapeDeriv[0][2] = x42*y32-x32*y42;
  shapeDeriv[1][2] = x31*y43-x34*y13;
  shapeDeriv[2][2] = x24*y14-x14*y24;
  shapeDeriv[3][2] = x13*y21-x12*y31;
}
void femTetra4::evalJacobianMatrix(double coord1, double coord2, double coord3, femDoubleMat shDerivs){

}
double femTetra4::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  // Define Quantities
  double x21 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x32 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x43 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[2]]->coords[0];

  double y12 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y23 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y34 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];

  double z12 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z23 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[2]]->coords[2];
  double z34 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[3]]->coords[2];

  return (x21 * (y23 * z34 - y34 * z23) + x32 * (y34 * z12 - y12 * z34) + x43 * (y12 * z23 - y23 * z12));
}

femIntegrationRule* femTetra4::evalIntegrationRule(intRuleType type){
  femIntegrationRule* res = new femIntegrationRule(type);
  return res;
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
void femTetra10::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femTetra10::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,
                                             femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}
void femTetra10::evalJacobianMatrix(double coord1, double coord2, double coord3, femDoubleMat shDerivs){
  throw femException("Not Implemented.\n");
}
double femTetra10::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  throw femException("Not Implemented.\n");
}
femIntegrationRule* femTetra10::evalIntegrationRule(intRuleType type){
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
double femTetra4::EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){
  double minMixedProduct = std::numeric_limits<double>::max();
  double vecA[3] = {0.0};
  double vecB[3] = {0.0};
  double vecC[3] = {0.0};
  double vec1[3] = {0.0};
  int n0,n1,n2,n3;
  double currProd = 0.0;
  for(unsigned int loopA=0;loopA<elementConnections.size();loopA++){
    // Get Local Nodes
    n0 = elementConnections[(loopA) % (elementConnections.size())];
    n1 = elementConnections[(loopA+1) % (elementConnections.size())];
    n2 = elementConnections[(loopA+2) % (elementConnections.size())];
    n3 = elementConnections[(loopA+3) % (elementConnections.size())];
    // Get three vectors
    vecA[0] = nodeList[n1]->coords[0] - nodeList[n0]->coords[0];
    vecA[1] = nodeList[n1]->coords[1] - nodeList[n0]->coords[1];
    vecA[2] = nodeList[n1]->coords[2] - nodeList[n0]->coords[2];
    vecB[0] = nodeList[n2]->coords[0] - nodeList[n0]->coords[0];
    vecB[1] = nodeList[n2]->coords[1] - nodeList[n0]->coords[1];
    vecB[2] = nodeList[n2]->coords[2] - nodeList[n0]->coords[2];
    vecC[0] = nodeList[n3]->coords[0] - nodeList[n0]->coords[0];
    vecC[1] = nodeList[n3]->coords[1] - nodeList[n0]->coords[1];
    vecC[2] = nodeList[n3]->coords[2] - nodeList[n0]->coords[2];
    // Get external product
    femUtils::Do3DExternalProduct(vecA,vecB,vec1);
    // Get Internal Product
    currProd = femUtils::Do3DInternalProduct(vec1,vecC);
    // Store Minimum Value
    if(currProd < minMixedProduct){
      minMixedProduct = currProd;
    }
  }
  return minMixedProduct;
}

// TO IMPLEMENT
double femTetra10::EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){
  throw femException("Internal: femTetra10::EvalMixProduct not Implemented");
  return 0.0;
}

// TO IMPLEMENT
double femTri3::EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){
  throw femException("Internal: femTria3::EvalMixProduct not Implemented");
  return 0.0;
}


// ================================
// FINITE ELEMENT ROUTINES FOR TRI3
// ================================
// =====================================
// VIRTUAL FUNCTIONS FOR FINITE ELEMENTS
// =====================================
void femTri3::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  throw femException("Not Implemented.\n");
}
void femTri3::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleMat &shapeDeriv){
  throw femException("Not Implemented.\n");
}
void femTri3::evalJacobianMatrix(double coord1, double coord2, double coord3, femDoubleMat shDerivs){
  throw femException("Not Implemented.\n");
}
double femTri3::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  throw femException("Not Implemented.\n");
}
femIntegrationRule* femTri3::evalIntegrationRule(intRuleType type){
  throw femException("Not Implemented.\n");
}

// ====================================================
// ASSEMBLE MASS MATRIX FOR ADVECTION DIFFUSION PROBLEM
// ====================================================
void femTetra4::assembleMass(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double** massMat){
  femDoubleVec shapeFunction;
  femDoubleMat shapeFunctionDerivs;
  int globNode = 0;
  double currProd = 0.0;
  for(int loopA=0;loopA<rule.totalGP;loopA++){
    // Eval Shape Functions
    evalShapeFunction(nodeList,rule.coords[loopA][0],rule.coords[loopA][1],rule.coords[loopA][2],shapeFunction);
    // Eval Shape Function Derivatives
    evalShapeFunctionDerivative(nodeList,rule.coords[loopA][0],rule.coords[loopA][1],rule.coords[loopA][2],shapeFunctionDerivs);
    // Loop Through the Nodes
    for(int loopB=0;loopB<numberOfNodes;loopB++){
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        // Eval Product with Velocity
        currProd = 0.0;
        for(int loopD=0;loopD<kDims;loopD++){
          globNode = elementConnections[loopB];
          currProd += nodeVelocities[globNode][loopD] * shapeFunctionDerivs[loopB][loopD];
        }
        massMat[loopB][loopC] += shapeFunction[loopB]*shapeFunction[loopC] + tauSUPG[loopA] * currProd * shapeFunction[loopC] * rule.weight[loopA];
      }
    }
  }
}

// =========================================================
// ASSEMBLE STIFFNESS MATRIX FOR ADVECTION DIFFUSION PROBLEM
// =========================================================
void femTetra4::assembleStiffness(femDoubleMat &nodeVelocities, std::vector<femNode*> nodeList, std::vector<double> tauSUPG, femIntegrationRule rule, double diffusivity, femDoubleMat &stiffnessMat){
  femDoubleVec shapeFunction;
  femDoubleMat shapeFunctionDerivs;
  int globNode = 0;
  double uNablaNB = 0.0;
  double uNablaNA = 0.0;
  double NablaNANablaNB = 0.0;
  double firstTerm = 0.0;
  double secondTerm = 0.0;
  double tauTerm = 0.0;
  for(int loopA=0;loopA<rule.totalGP;loopA++){
    // Eval Shape Functions
    evalShapeFunction(nodeList,rule.coords[loopA][0],rule.coords[loopA][1],rule.coords[loopA][2],shapeFunction);
    // Eval Shape Function Derivatives
    evalShapeFunctionDerivative(nodeList,rule.coords[loopA][0],rule.coords[loopA][1],rule.coords[loopA][2],shapeFunctionDerivs);
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
        stiffnessMat[loopB][loopC] += firstTerm + secondTerm + tauTerm;
      }
    }
  }
}
