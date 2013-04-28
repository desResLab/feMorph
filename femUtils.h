#ifndef FEMUTILS_H
#define FEMUTILS_H

#include "femConstants.h"
#include "femNode.h"
#include "femInputData.h"
#include "femModelSlice.h"

#include <vector>
#include <math.h>

// Random Number Generator
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace femUtils{

// GLOBAL VARIABLE FOR GENERATION
static boost::random::mt19937 intGen;

// Generate Uniform Integer (global Variable Defined)
inline int GenerateUniformIntegers(int lowIdx, int upIdx) {
    boost::random::uniform_int_distribution<> dist(lowIdx, upIdx);
    return dist(intGen);
}

// =============
// Write Message
// =============
inline void WriteMessage(std::string m){
   printf("%s",m.c_str());
   fflush(stdout);
}

// ========================
// Write Application Header
// ========================
inline void WriteAppHeader(){
  WriteMessage("\n");
  WriteMessage("------------------------------------\n");
  WriteMessage("FEM Morphing Application\n");
  WriteMessage("Release 0.1 beta\n");
  WriteMessage("Daniele Schiavazzi Ph.D., 2013\n");
  WriteMessage("Marsden LAB\n");
  WriteMessage("University of California @ San Diego\n");
  WriteMessage("------------------------------------\n");
  WriteMessage("\n");
}

// ======================
// Write Application Help
// ======================
inline void WriteAppHelp(){
  WriteMessage("Usage: feMorph [OPTION] [FILE]\n");
  WriteMessage("Mandatory option arguments:\n");
  WriteMessage("-n               Normal Execution\n");
  WriteMessage("-d               Debug Mode, a number of model and reference entities are \n");
  WriteMessage("                 printer using the VTK legacy file format. \n");
  WriteMessage("-?,-h            this help\n");
}

// Extral product between vectors
inline void Do3DExternalProduct(double* v1,double* v2,double* v3){
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v2[2]*v1[0];
  v3[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

// Eval Internal Product
inline double Do3DInternalProduct(double* v1, double* v2){
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// ===================
// Normalize 3D Vector
// ===================
inline void Normalize3DVector(double* vector){
  // Define modulus
  double modulus = sqrt((vector[0]*vector[0])+(vector[1]*vector[1])+(vector[2]*vector[2]));
  // Normalize
  if (modulus>0.0){
    vector[0] = vector[0]/modulus;
  }
  if (modulus>0.0){
    vector[1] = vector[1]/modulus;
  }
  if (modulus>0.0) {
    vector[2] = vector[2]/modulus;
  }
}

// =================================================
// General Rotation Around a Vector with quaternions
// The Angle is in degrees
// =================================================
inline void Rotate3DVectorAroundAxis(double* vec, double angle,double* axis){
  double quatMat[4][4];
  double quatVect[4] = {0.0};
  double quatRes[4] = {0.0};
  // Normalize Axis
  Normalize3DVector(axis);
  // Write the compoents of the quaternion
  double w = cos(0.5*angle*(kPI/double(180.0)));
  double x = axis[0]*sin(0.5*angle*(kPI/double(180.0)));
  double y = axis[1]*sin(0.5*angle*(kPI/double(180.0)));
  double z = axis[2]*sin(0.5*angle*(kPI/double(180.0)));
  // Eval rotation Matrix
  quatMat[0][0] = w*w+x*x-y*y-z*z;
  quatMat[0][1] = 2.0*x*y+2.0*w*z;
  quatMat[0][2] = 2.0*x*z-2.0*w*y;
  quatMat[0][3] = 0.0;
  quatMat[1][0] = 2.0*x*y-2.0*w*z;
  quatMat[1][1] = w*w-x*x+y*y-z*z;
  quatMat[1][2] = 2.0*y*z+2.0*w*x;
  quatMat[1][3] = 0.0;
  quatMat[2][0] = 2.0*x*z+2.0*w*y;
  quatMat[2][1] = 2.0*y*z-2.0*w*x;
  quatMat[2][2] = w*w-x*x-y*y+z*z;
  quatMat[2][3] = 0.0;
  quatMat[3][0] = 0.0;
  quatMat[3][1] = 0.0;
  quatMat[3][2] = 0.0;
  quatMat[3][3] = w*w+x*x+y*y+z*z;
  // Expand vector to rotate
  quatVect[3] = 0.0;
  for(int loopA=0;loopA<3;loopA++){
    quatVect[loopA] = vec[loopA];
  }
  // Matrix vector multiplication
  for(int loopA=0;loopA<4;loopA++){
    quatRes[loopA] = 0.0;
    for(int loopB=0;loopB<4;loopB++){
      quatRes[loopA] += quatVect[loopB] * quatMat[loopB][loopA];
    }
  }
  // Final Copy of the Result
  for(int loopA=0;loopA<3;loopA++){
    vec[loopA] = quatRes[loopA];
  }
}

// ============================================
// Export Stenosys Rotate and Scaled Box to VTK
// ============================================
inline void ExportStenosisBoxToVTK(std::string fileName, std::vector<femNode*> steNodeList){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");

  // Wrtie Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Model Exported from feMorph\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");

  // Write Points
  fprintf(outFile,"POINTS %d float\n",int(steNodeList.size()));
  for(unsigned int loopA=0;loopA<steNodeList.size();loopA++){
    fprintf(outFile,"%e %e %e\n",steNodeList[loopA]->coords[0],steNodeList[loopA]->coords[1],steNodeList[loopA]->coords[2]);
  }

  // Write single CELL header
  fprintf(outFile,"CELLS %d %d\n",1,9);
  fprintf(outFile,"8 0 1 2 3 4 5 6 7\n");

  // Write Cells Type Header
  fprintf(outFile,"CELL_TYPES 1\n");
  fprintf(outFile,"11\n");

  // Close Output file
  fclose(outFile);
}

// =============================================
// Export Reference Axis System to VTK for Debug
// =============================================
inline void ExportReferenceToVTKLegacy(femInputData* data, std::string fileName){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");

  // Write Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Model Exported from feMorph\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");

  double currDisps[6] = {0.0};
  double currCoord[3] = {0.0};
  std::vector<femNode*> nodeList;

  // Include the origin
  currCoord[0] = data->mainModelOrigin[0];
  currCoord[1] = data->mainModelOrigin[1];
  currCoord[2] = data->mainModelOrigin[2];
  femNode* newNode = new femNode(0,currCoord,currDisps);
  nodeList.push_back(newNode);

  // Loop through the three main axis
  double axis[3] = {0.0};
  for(int loopA=0;loopA<3;loopA++){
    axis[0] = data->mainModelRefSystem[0][loopA];
    axis[1] = data->mainModelRefSystem[1][loopA];
    axis[2] = data->mainModelRefSystem[2][loopA];
    Normalize3DVector(axis);
    // Positive
    currCoord[0] = data->mainModelOrigin[0] + 0.5*data->stenosisLength*axis[0];
    currCoord[1] = data->mainModelOrigin[1] + 0.5*data->stenosisLength*axis[1];
    currCoord[2] = data->mainModelOrigin[2] + 0.5*data->stenosisLength*axis[2];
    newNode = new femNode(loopA*2+1,currCoord,currDisps);
    nodeList.push_back(newNode);
    // Negative
    currCoord[0] = data->mainModelOrigin[0] - 0.5*data->stenosisLength*axis[0];
    currCoord[1] = data->mainModelOrigin[1] - 0.5*data->stenosisLength*axis[1];
    currCoord[2] = data->mainModelOrigin[2] - 0.5*data->stenosisLength*axis[2];
    newNode = new femNode(loopA*2+2,currCoord,currDisps);
    nodeList.push_back(newNode);
  }

  // Write Points
  fprintf(outFile,"POINTS %d float\n",int(nodeList.size()));
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    fprintf(outFile,"%e %e %e\n",nodeList[loopA]->coords[0],nodeList[loopA]->coords[1],nodeList[loopA]->coords[2]);
  }

  // Write single CELL header
  fprintf(outFile,"CELLS 6 18\n");
  fprintf(outFile,"2 0 1\n");
  fprintf(outFile,"2 0 2\n");
  fprintf(outFile,"2 0 3\n");
  fprintf(outFile,"2 0 4\n");
  fprintf(outFile,"2 0 5\n");
  fprintf(outFile,"2 0 6\n");

  // Write Cells Type Header
  fprintf(outFile,"CELL_TYPES 6\n");
  fprintf(outFile,"4\n");
  fprintf(outFile,"4\n");
  fprintf(outFile,"4\n");
  fprintf(outFile,"4\n");
  fprintf(outFile,"4\n");
  fprintf(outFile,"4\n");

  // Close Output file
  fclose(outFile);
}

// Get Average Nodee From List
inline void GetAverageNodeFromList(std::vector<femNode*> &nodeList,double* stenosisBoxCenter){
  stenosisBoxCenter[0] = 0.0;
  stenosisBoxCenter[1] = 0.0;
  stenosisBoxCenter[2] = 0.0;
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    for(int loopB=0;loopB<3;loopB++){
      stenosisBoxCenter[loopB] += nodeList[loopA]->coords[loopB];
    }
  }
  // Make Average
  for(int loopB=0;loopB<3;loopB++){
    stenosisBoxCenter[loopB] /= double(nodeList.size());
  }
}

// symmetric round down
// Bias: towards zero
template <typename FloatType>
inline FloatType trunc(const FloatType& value)
{
  FloatType result = std::floor(std::fabs(value));
  return (value < 0.0) ? -result : result;
}

// =================================
// Do the Euclidean Norm of a Vector
// =================================
inline double DoEucNorm(double* vec){
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

// ====================================================
// Find intersection between coordinate plane and point
// ====================================================
inline void findPointsToCoordinate1PlaneIntersection(double planeCoord, double* coords1, double* coords2, std::vector<femNode*> &intNodeList){
  double param = 0.0;
  double intCoord[3] = {0.0};
  femNode* node;
  if(fabs(coords2[0] - coords1[0])>kMathZero){
    param = (planeCoord - coords1[0])/(coords2[0] - coords1[0]);
    if((param>=0.0)&&(param<=1.0)){
      // Get Intersection
      intCoord[0] = coords1[0] + param*(coords2[0] - coords1[0]);
      intCoord[1] = coords1[1] + param*(coords2[1] - coords1[1]);
      intCoord[2] = coords1[2] + param*(coords2[2] - coords1[2]);
      node = new femNode(0,intCoord[0],intCoord[1],intCoord[2]);
      intNodeList.push_back(node);
    }
  }else if(fabs(planeCoord-coords2[0])<kMathZero){
    // The side is on the plane: insert both points
    // Point 1
    node = new femNode(0,intCoord[0],intCoord[1],intCoord[2]);
    intNodeList.push_back(node);
    // Point 2
    node = new femNode(0,intCoord[0],intCoord[1],intCoord[2]);
    intNodeList.push_back(node);
  }
}


// ================
// Eval Intersetion
// ================
inline void evalCoordinatePlaneIntersections(double currSliceCoord, double *coords1, double* coords2, double* coords3, std::vector<femNode*> &intNodeList){
  // 1 - 2
  findPointsToCoordinate1PlaneIntersection(currSliceCoord,coords1,coords2,intNodeList);
  // 2 - 3
  findPointsToCoordinate1PlaneIntersection(currSliceCoord,coords2,coords3,intNodeList);
  // 3 - 1
  findPointsToCoordinate1PlaneIntersection(currSliceCoord,coords3,coords1,intNodeList);
}

// =========================
// Plot Slices to VTK legacy
// =========================
inline void PlotSlicesToVTK(std::string fileName, std::vector<femModelSlice*> &slices){

  // Write Message
  femUtils::WriteMessage(std::string("(debug) Exporting Model Slices to VTK..."));

  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");

  // Write Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"feMorph Slice Data\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");

  // sum up point total
  int pointTotal = 0;
  for(unsigned int loopA=0;loopA<slices.size();loopA++){
    pointTotal += slices[loopA]->slicedPoints.size();
  }

  // Write Points
  fprintf(outFile,"POINTS %d float\n",pointTotal);
  for(unsigned int loopA=0;loopA<slices.size();loopA++){
    for(unsigned int loopB=0;loopB<slices[loopA]->slicedPoints.size();loopB++){
      fprintf(outFile,"%e %e %e\n",slices[loopA]->slicedPoints[loopB]->coords[0],
                                   slices[loopA]->slicedPoints[loopB]->coords[1],
                                   slices[loopA]->slicedPoints[loopB]->coords[2]);
    }
  }

  // Write CELLS header
  int firstPointNumber = 0;
  int currentPointNumber = 0;
  fprintf(outFile,"CELLS %d %d\n",pointTotal,3*pointTotal);
  for(unsigned int loopA=0;loopA<slices.size();loopA++){
    for(unsigned int loopB=0;loopB<slices[loopA]->slicedPoints.size();loopB++){
      if(loopB<slices[loopA]->slicedPoints.size()-1){
        // Before the last point
        fprintf(outFile,"2 %d %d\n",currentPointNumber,currentPointNumber+1);
      }else{
        // Last point
        fprintf(outFile,"2 %d %d\n",currentPointNumber,firstPointNumber);

      }
      // Increment Current Point Number
      currentPointNumber++;
    }
    // Update Counter
    firstPointNumber += slices[loopA]->slicedPoints.size();
  }

  // Write Cells Type Header
  fprintf(outFile,"CELL_TYPES %d\n",pointTotal);
  for(int loopA=0;loopA<pointTotal;loopA++){
    fprintf(outFile,"4\n");
  }

  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));

  // Close Output file
  fclose(outFile);
}


}

#endif // FEMUTILS_H

