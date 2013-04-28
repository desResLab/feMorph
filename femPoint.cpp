#include "femPoint.h"

// Constructor
femPoint::femPoint(double* localCoords){
  // Assign Coords
  coords[0] = localCoords[0];
  coords[1] = localCoords[1];
  coords[2] = localCoords[2];
}

// Distructor
femPoint::~femPoint(){

}

