#include "femPointGrid.h"

// Create a regular grid
femPointGrid::femPointGrid(double* lim, int* totals){
  double spacing[3];
  // Assign Limits
  for(int loopA=0;loopA<6;loopA++){
    limits[loopA] = lim[loopA];
  }
  // Assign Totals
  totPoints[0] = totals[0];
  totPoints[1] = totals[1];
  totPoints[2] = totals[2];
  // Assign Spacing
  spacing[0] = limits[2*0+1] - limits[2*0+1]/(double)totPoints[0];
  spacing[1] = limits[2*1+1] - limits[2*1+1]/(double)totPoints[1];
  spacing[2] = limits[2*2+1] - limits[2*2+1]/(double)totPoints[2];

  int totalPoints = totPoints[0] * totPoints[1] * totPoints[2];

  coords.resize(totalPoints);
  displacements.resize(totalPoints);
  for(int loopA=0;loopA<totalPoints;loopA++){
    coords[loopA].resize(3);
    displacements[loopA].resize(3);
    coords[loopA][0] = 0.0;
    coords[loopA][1] = 0.0;
    coords[loopA][2] = 0.0;
    displacements[loopA][0] = 0.0;
    displacements[loopA][1] = 0.0;
    displacements[loopA][2] = 0.0;
  }

  int count = 0;
  for(int loopA=0;loopA<totPoints[2];loopA++){
    for(int loopB=0;loopB<totPoints[1];loopB++){
      for(int loopC=0;loopC<totPoints[0];loopC++){
        coords[count][0] = limits[2*0] + loopC * spacing[0];
        coords[count][1] = limits[2*1] + loopB * spacing[1];
        coords[count][2] = limits[2*2] + loopA * spacing[2];
        count++;
      }
    }
  }
}

femPointGrid::femPointGrid(double* limits, int* totals, const femDoubleMat& initCoords){
  // Call Standard constructor
  femPointGrid(limits,totals);
  // Replace Coords
  for(int loopA=0;loopA<totalPoints;loopA++){
    // Get the Coords of the current point
    coords[loopA][0] = initCoords[loopA][0];
    coords[loopA][1] = initCoords[loopA][1];
    coords[loopA][2] = initCoords[loopA][2];
  }
}


// Member Functions
// Evaluate the deformations at a set of point locations
femPointGrid::morphPoints(const femDoubleMat& pointLocations, femDoubleMat& pointDeformations){

}

