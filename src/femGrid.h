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
  virtual ~femGrid();
  // Member Functions
  // To Single Index
  int ToIndexes(double* coord);
  // To Index Array
  void ToIndexArray(double* coord,int* index);
  // From inde Array to single index
  int IndexArrayToIndex(int idx0,int idx1,int idx2);
  void ExportToVTKLegacy(std::string fileName);
};

#endif // FEMGRID_H
