#include "femFace.h"

femFace::femFace()
{
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
