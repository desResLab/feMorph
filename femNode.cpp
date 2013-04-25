#include "femNode.h"

// Constructor for femNode
femNode::femNode(int number, double coordX, double coordY, double coordZ){
  // Assign Node Number
  nodeNumber = number;
  // Assign Node Coordinates
  coords[0] = coordX;
  coords[1] = coordY;
  coords[2] = coordZ;
  // Initialize Displacements
  for(int loopA=0;loopA<6;loopA++){
    displacements[loopA] = 0.0;
  }
}

// Other Constructor
femNode::femNode(int number, double* coords, double* disps){
  // Assign Node Number
  nodeNumber = number;
  // Assign Node Coordinates
  for(int loopA=0;loopA<3;loopA++){
    this->coords[loopA] = coords[loopA];
  }
  for(int loopA=0;loopA<6;loopA++){
    this->displacements[loopA] = disps[loopA];
  }
}

// Set Node Displacements
void femNode::setDisplacements(double DX, double DY, double DZ, double RX, double RY, double RZ){
  // Assign Displacements
  displacements[0] = DX;
  displacements[1] = DY;
  displacements[2] = DZ;
  displacements[3] = RX;
  displacements[4] = RY;
  displacements[5] = RZ;
}

// ===========================================
// Transform Nodes coordinates to a new System
// ===========================================
void femNode::TransformNodeCoords(double* origin, double** rotMat, double* newCoords){
  // Allocate
  double tempCoords[3] = {0.0};
  // Set to the intial Coords plus translation
  for(int loopA=0;loopA<3;loopA++){
    tempCoords[loopA] = coords[loopA] - origin[loopA];
  }
  // Rotate
  for(int loopA=0;loopA<3;loopA++){
    newCoords[loopA] = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      newCoords[loopA] += rotMat[loopB][loopA] * tempCoords[loopB];
    }
    // Add origin: NO!
    // newCoords[loopA] += origin[loopA];
  }
}

// =============================================
// Transform Nodes displacements to a new System
// =============================================
void femNode::TransformDisplacements(double** rotMat, double* newDisps){
  // Allocate
  double tempCoords[3] = {0.0};
  // Set to the intial Coords plus translation
  for(int loopA=0;loopA<3;loopA++){
    tempCoords[loopA] = coords[loopA];
  }
  // Rotate
  for(int loopA=0;loopA<3;loopA++){
    newDisps[loopA] = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      newDisps[loopA] += rotMat[loopA][loopB] * tempCoords[loopB];
    }
  }
}
