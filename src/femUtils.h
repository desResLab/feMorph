#ifndef FEMUTILS_H
#define FEMUTILS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <vector>

#include "femConstants.h"
#include "femNode.h"
#include "femInputData.h"
#include "femModelSlice.h"
#include "femException.h"
#include "femWeightedFileName.h"
#include "femTypes.h"

// Random Number Generator
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/algorithm/string.hpp>

namespace femUtils{

// RANDOM NUMBERS
static boost::random::mt19937 intGen;
int GenerateUniformIntegers(int lowIdx, int upIdx);

// APP MESSAGES 
void WriteMessage(std::string m);
void WriteAppHeader();
void WriteAppHelp();

// ALBEGRA AND VECTOR/MATRIX OPERATIONS
void Do3DExternalProduct(double* v1,double* v2,double* v3);
double Do3DInternalProduct(double* v1, double* v2);
double Do3DMixedProduct(double* v1, double* v2, double* v3);
void Normalize3DVector(double* vector);
void matZeros(femDoubleMat& mat,int row_size,int column_size);
void Rotate3DVectorAroundAxis(double* vec, double angle,double* axis);
double DoEucNorm(double* vec);
double DoEucNorm(const femDoubleVec& vec);
void invert3x3Matrix(const femDoubleMat& m,femDoubleMat &invMat, double& detJ);
void invert3x3MatrixFor1DElements(const femDoubleMat& mat,femDoubleMat& invMat, double& detJ);
double getMaxModule(const femDoubleVec& vec);
double getMaxModule(const femDoubleMat& mat, int dim1, int dim2);
double coth(double value);
bool isInsideLimits(double* centroid,double* limitBox);
double getMatrixNorm(const femDoubleMat& mat);

// LIST INDICES AND SORTING
int findVectorInverseMapping(int value, const std::vector<int>& list);
void bubbleSortIntVector(std::vector<int>& nodes);
void MakeCompactList(const std::vector<int>& list, std::vector<int> &compactList);
void MakeInverseList(std::vector<int> list,int totalNodes,std::vector<int> &inverseList);

// TRUNCATION
int trunc(double value);

// VTK-RELATED OPERATIONS
void ExportStenosisBoxToVTK(std::string fileName, std::vector<femNode*> steNodeList);
void ExportReferenceToVTKLegacy(femInputData* data, std::string fileName);
void PlotSlicesToVTK(std::string fileName, std::vector<femModelSlice*> &slices);

// OTHER GEOMETRIC UTILITIES
void GetAverageNodeFromList(std::vector<femNode*> &nodeList,double* stenosisBoxCenter);
void findPointsToCoordinate1PlaneIntersection(double planeCoord, double* coords1, double* coords2, std::vector<femNode*> &intNodeList);
void evalCoordinatePlaneIntersections(double currSliceCoord, double *coords1, double* coords2, double* coords3, std::vector<femNode*> &intNodeList);
bool AreCompatibleBoxes(double* firstBox, double* secondBox,double tolerance);

// STRING OPERATIONS
string intToStr(int a);
string floatToStr(float a);
std::string extractFileName(std::string fullPath);
bool file_exists(const string& fileName);

// INTERPOLATION
femDoubleVec interpVelocity(const double time, const double time1, const femDoubleVec& vel1, const double time2, const femDoubleVec& vel2,  int dims=3);
femDoubleVec interpTableData(double time, const femDoubleMat& table, int time_idx=3);

// I/O OPERATIONS
void ReadListFromFile(std::string fileName, std::vector<std::string> &fileList1);
void WriteGraphToFile(std::string fileName, int vecSize, std::vector<double> &vecX, std::vector<double> &vecY);
void writeVectorToFile(std::string fileName, const femDoubleVec& vec);
void writeMatrixToFile(std::string fileName, const femDoubleMat& mat);
void printMatrix(const femDoubleMat& mat);
void printVector(const femDoubleVec& vec);
void ReadWeightedFileList(std::string fileName, std::vector<femWeightedFileName*> &fileList);
}

#endif // FEMUTILS_H

