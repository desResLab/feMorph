#include "femModelSlice.h"
#include "femUtils.h"

femModelSlice::femModelSlice()
{
}

femModelSlice::~femModelSlice()
{
}

// ===============================================
// Sort angle and at the same time the permutation
// ===============================================
void sortAnglesAndPerm(std::vector<double> angles, std::vector<femPoint*> &slicedPoints){
  double tempAngle = 0.0;
  double tempCoords[3] = {0.0};
  for(unsigned int loopA=0;loopA<angles.size()-1;loopA++){
    for(unsigned int loopB=loopA;loopB<angles.size();loopB++){
      if(angles[loopA] > angles[loopB]){
        // Swap angles and permutations
        tempAngle = angles[loopB];
        angles[loopB] = angles[loopA];
        angles[loopA] = tempAngle;
        // Swap Points
        tempCoords[0] = slicedPoints[loopB]->coords[0];
        tempCoords[1] = slicedPoints[loopB]->coords[1];
        tempCoords[2] = slicedPoints[loopB]->coords[2];
        slicedPoints[loopB]->coords[0] = slicedPoints[loopA]->coords[0];
        slicedPoints[loopB]->coords[1] = slicedPoints[loopA]->coords[1];
        slicedPoints[loopB]->coords[2] = slicedPoints[loopA]->coords[2];
        slicedPoints[loopA]->coords[0] = tempCoords[0];
        slicedPoints[loopA]->coords[1] = tempCoords[1];
        slicedPoints[loopA]->coords[2] = tempCoords[2];
      }
    }
  }
}

// ================================
// Get Centroid of the slice points
// ================================
void femModelSlice::getPointCentoid(double* centerPoint){
  // Init
  centerPoint[0] = 0.0;
  centerPoint[1] = 0.0;
  centerPoint[2] = 0.0;
  for(unsigned int loopA=0;loopA<slicedPoints.size();loopA++){
    centerPoint[0] += slicedPoints[loopA]->coords[0];
    centerPoint[1] += slicedPoints[loopA]->coords[1];
    centerPoint[2] += slicedPoints[loopA]->coords[2];
  }
  centerPoint[0] /= (double)slicedPoints.size();
  centerPoint[1] /= (double)slicedPoints.size();
  centerPoint[2] /= (double)slicedPoints.size();
}

// ===============
// Eval Slice Area
// ===============
double femModelSlice::EvalSliceArea(double** refSystem){
  double axisSign[3];
  double refModelAxis[3];
  double currentSign = 0.0;
  refModelAxis[0] = refSystem[0][0];
  refModelAxis[1] = refSystem[1][0];
  refModelAxis[2] = refSystem[2][0];
  double centerPoint[3] = {0.0};
  double refVec[3] = {0.0};
  double otherVec[3]  = {0.0};
  double vec1[3] = {0.0};
  double vec2[3] = {0.0};
  double vec3[3] = {0.0};
  // Get Point centroid
  getPointCentoid(centerPoint);
  // Build Reference Vector
  refVec[0] = centerPoint[0] - slicedPoints[0]->coords[0];
  refVec[1] = centerPoint[1] - slicedPoints[0]->coords[1];
  refVec[2] = centerPoint[2] - slicedPoints[0]->coords[2];
  femUtils::Normalize3DVector(refVec);
  std::vector<double> angles;
  angles.push_back(0.0);
  // Eval the angle with reference Vector
  double currAngle = 0.0;
  for(unsigned int loopA=1;loopA<slicedPoints.size();loopA++){
    otherVec[0] = centerPoint[0] - slicedPoints[loopA]->coords[0];
    otherVec[1] = centerPoint[1] - slicedPoints[loopA]->coords[1];
    otherVec[2] = centerPoint[2] - slicedPoints[loopA]->coords[2];
    // Get Angle
    femUtils::Normalize3DVector(otherVec);
    currAngle = acos(femUtils::Do3DInternalProduct(refVec,otherVec))*(180.0/kPI);
    // Get Sign
    femUtils::Do3DExternalProduct(refVec,otherVec,axisSign);
    femUtils::Normalize3DVector(axisSign);
    currentSign = femUtils::Do3DInternalProduct(refModelAxis,axisSign);
    angles.push_back(currentSign*currAngle);
  }
  // Sort
  sortAnglesAndPerm(angles,slicedPoints);
  // Scan triangle area
  double area = 0.0;
  double l1 = 0.0;
  double l2 = 0.0;
  double l3 = 0.0;
  double s = 0.0;
  int firstIndex = 0;
  int secondIndex = 0;
  for(unsigned int loopA=0;loopA<slicedPoints.size();loopA++){
    if(loopA == (slicedPoints.size()-1)){
      firstIndex = slicedPoints.size()-1;
      secondIndex = 0;
    }else{
      firstIndex = loopA;
      secondIndex = (loopA + 1);
    }
    vec1[0] = centerPoint[0] - slicedPoints[firstIndex]->coords[0];
    vec1[1] = centerPoint[1] - slicedPoints[firstIndex]->coords[1];
    vec1[2] = centerPoint[2] - slicedPoints[firstIndex]->coords[2];
    l1 = femUtils::DoEucNorm(vec1);
    vec2[0] = centerPoint[0] - slicedPoints[secondIndex]->coords[0];
    vec2[1] = centerPoint[1] - slicedPoints[secondIndex]->coords[1];
    vec2[2] = centerPoint[2] - slicedPoints[secondIndex]->coords[2];
    l2 = femUtils::DoEucNorm(vec2);
    vec3[0] = slicedPoints[firstIndex]->coords[0] - slicedPoints[secondIndex]->coords[0];
    vec3[1] = slicedPoints[firstIndex]->coords[1] - slicedPoints[secondIndex]->coords[1];
    vec3[2] = slicedPoints[firstIndex]->coords[2] - slicedPoints[secondIndex]->coords[2];
    l3 = femUtils::DoEucNorm(vec3);
    s = 0.5*(l1+l2+l3);
    // sum areas up with Heron's formula
    area += sqrt(s*(s-l1)*(s-l2)*(s-l3));
  }
  // Return
  return area;
}
