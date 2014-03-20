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
femElement::~femElement()
{
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
  //  throw new femException("Internal: Tet4 Volume Coordinates do not sum up to one.");
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
  //  throw new femException("Error: Internal. Tet4 Shape Function don't sum up to one");
  //}
  // Return
  return isInside;
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
  //  throw new femException("Error: Internal. In Tet10 internal Check, Shape Function don't sum up to one");
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
      throw new femException("Internal: Invalid Element connections in face.\n");
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
  throw new femException("Internal: femTetra10::EvalMixProduct not Implemented");
  return 0.0;
}

// TO IMPLEMENT
double femTri3::EvalMixProduct(double dispFactor, std::vector<femNode*> &nodeList){
  throw new femException("Internal: femTria3::EvalMixProduct not Implemented");
  return 0.0;
}
