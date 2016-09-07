#include "femFace.h"
#include "femUtils.h"

femFace::femFace()
{
  // Initialize Group to -1
  group = -1;
}

// Aternative Constructor With Node List
femFace::femFace(std::vector<int> nodes){
  for(unsigned int loopA=0;loopA<nodes.size();loopA++){
    faceNodes.push_back(nodes[loopA]);
  }
}

// Alternative Constructor
femFace::femFace(int tempNumber, std::vector<int> nodes){
  number = tempNumber;
  for(unsigned int loopA=0;loopA<nodes.size();loopA++){
    faceNodes.push_back(nodes[loopA]);
  }
}

femFace::femFace(femFace* other){
  // Face Nodes
  for(unsigned int loopA=0;loopA<other->faceNodes.size();loopA++){
    faceNodes.push_back(other->faceNodes[loopA]);
  }
  // Face Elements
  for(unsigned int loopA=0;loopA<other->faceElements.size();loopA++){
    faceElements.push_back(other->faceElements[loopA]);
  }
}

// ==================
// EVAL FACE CENTROID
// ==================
void femFace::evalFaceCentroid(std::vector<femNode*> nodeList, double* centroid){
  int currNode = 0;
  for(int loopA=0;loopA<kDims;loopA++){
    centroid[loopA] = 0.0;
  }
  for(size_t loopA=0;loopA<faceNodes.size();loopA++){
    currNode = faceNodes[loopA];
    for(int loopA=0;loopA<kDims;loopA++){
      centroid[loopA] = centroid[loopA] + nodeList[currNode]->coords[loopA];
    }
  }
  // Eval Average
  for(int loopA=0;loopA<kDims;loopA++){
    centroid[loopA] /= (double)faceNodes.size();
  }
}

// ================
// EVAL FACE NORMAL
// ================
void femFace::evalFaceNormal(std::vector<femNode*> nodeList, double* normal){
  int node1 = faceNodes[0];
  int node2 = faceNodes[1];
  int node3 = faceNodes[2];
  double firstVec[3] = {0.0};
  double secondVec[3] = {0.0};
  for(int loopA=0;loopA<kDims;loopA++){
    firstVec[loopA] = nodeList[node2]->coords[loopA] - nodeList[node1]->coords[loopA];
    secondVec[loopA] = nodeList[node3]->coords[loopA] - nodeList[node1]->coords[loopA];
  }
  femUtils::Do3DExternalProduct(firstVec,secondVec,normal);
  femUtils::Normalize3DVector(normal);
}
