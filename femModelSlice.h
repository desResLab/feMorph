#ifndef FEMMODELSLICE_H
#define FEMMODELSLICE_H

#include <stdlib.h>
#include <vector>

#include "femPoint.h"

class femModelSlice
{
public:
  // Data members
  std::vector<int> slicedFaces;
  std::vector<femPoint*> slicedPoints;
  // Constructor and Distructor
  femModelSlice();
  ~femModelSlice();
  // Data Members
  double EvalSliceArea(double** refSystem);
  void getPointCentoid(double* centerPoint);
};

#endif // FEMMODELSLICE_H
