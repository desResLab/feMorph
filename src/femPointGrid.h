#ifndef FEMPOINTGRID_H
#define FEMPOINTGRID_H

# include "femTypes.h"

class femPointGrid{
public:
  // Data Member
  double limits[6];
  double totPoints[3];
  femDoubleMat initialCoords;
  femDoubleMat displacements;

  // Constructor - with no initial coords
  femPointGrid();
  // Constructor - with provided initial coords
  femPointGrid(const femDoubleMat& initCoords);
  // Distructor
  virtual ~femPointGrid();

  // Member Functions
  // Evaluate the deformations at a set of point locations
  morphPoints(const femDoubleMat& pointLocations, femDoubleMat& pointDeformations);
};

#endif // FEMPOINTGRID_H
