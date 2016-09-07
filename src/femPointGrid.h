#ifndef FEMPOINTGRID_H
#define FEMPOINTGRID_H

# include <math.h>
# include "femConstants.h"
# include "femUtils.h"
# include "femTypes.h"
# include "femModel.h"

class femPointGrid{
public:
  // Data Member
  double minPoint[3];
  double centerPoint[3];
  double size[3];
  int totPoints[3];
  double gridAxis[3][3];
  femDoubleMat coords;
  femDoubleMat disps;

  // Constructor - with no initial coords
  femPointGrid(double* locMinPoint,
               int* locTotPoints,
               double* locGridAxis_S,
               double* locGridAxis_T,
               double* locGridAxis_U,
               const femIntVec& dispNodes,
               const femDoubleMat& dispVals);
  // Distructor
  virtual ~femPointGrid();

  // Member Functions
  // Evaluate Trivariate Bernstein Polynomial
  void evalTrilinearBernsteinPolynomials(double* currCoord,double* currResult);
  // Evaluate the deformations at a set of point locations
  void morphPoints(const femDoubleMat& pointLocations, femDoubleMat& pointDeformations);
  // Transform Space Coords to Box Coords
  void getBoxCoords(const femDoubleMat& spaceCoords, femDoubleMat& boxCoords);
  // Morph Mesh
  void morphModel(femModel* model);
  // Export Grid to VTK Legacy
  void ExportToVTKLegacy(std::string fileName);


};

#endif // FEMPOINTGRID_H
