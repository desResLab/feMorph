#include "femConstants.h"
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
void femNode::TransformNodeCoords(double* origin, double** rotMat, double* newCoords, const int transformType, const int dispType, double dispFactor){
  // Allocate
  double temp0Coords[3] = {0.0};
  double temp1Coords[3] = {0.0};
  double temp2Coords[3] = {0.0};

  // Add displacements if required
  if((dispType == kDeformed)&&(transformType == kDirect)){
    for(int loopA=0;loopA<3;loopA++){
      temp0Coords[loopA] = coords[loopA] + dispFactor*displacements[loopA];
    }
  }else{
    for(int loopA=0;loopA<3;loopA++){
      temp0Coords[loopA] = coords[loopA];
    }
  }

  // Set to the intial Coords plus translation
  if(transformType == kDirect){
    for(int loopA=0;loopA<3;loopA++){
      temp1Coords[loopA] = temp0Coords[loopA] - origin[loopA];
    }
  }else{
    for(int loopA=0;loopA<3;loopA++){
      temp1Coords[loopA] = temp0Coords[loopA];
    }
  }

  // Rotate
  for(int loopA=0;loopA<3;loopA++){
    temp2Coords[loopA] = 0.0;
    for(int loopB=0;loopB<3;loopB++){
      if(transformType == kDirect){
        temp2Coords[loopA] += rotMat[loopB][loopA] * temp1Coords[loopB];
      }else{
        temp2Coords[loopA] += rotMat[loopA][loopB] * temp1Coords[loopB];
      }
    }
  }

  // Add origin if reverse
  if(transformType == kReverse){
    for(int loopA=0;loopA<3;loopA++){
      newCoords[loopA] = temp2Coords[loopA] + origin[loopA];
    }
  }else{
    for(int loopA=0;loopA<3;loopA++){
      newCoords[loopA] = temp2Coords[loopA];
    }
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
