#ifndef FEMGRID_H
#define FEMGRID_H

#include <stdio.h>
#include <vector>

#include "femModel.h"
#include "femGridCell.h"

class femGrid{
public:
  // Data Member
  double gridLimits[6];
  double gridSpacing[3];
  std::vector<femGridCell*> gridData;
  // Constructor and Distructor
  femGrid(femModel* model);
  ~femGrid();
  // Member Functions
  int ToIndexes(double* coord);
  void ExportToVTKLegacy(std::string fileName);
};

#endif // FEMGRID_H
