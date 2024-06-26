#ifndef FEMUTILS_CPP
#define FEMUTILS_CPP

#include "femUtils.h"

namespace femUtils{

// Generate Uniform Integer (global Variable Defined)
int GenerateUniformIntegers(int lowIdx, int upIdx) {
    boost::random::uniform_int_distribution<> dist(lowIdx, upIdx);
    return dist(intGen);
}

// =============
// Write Message
// =============
void WriteMessage(std::string m){
   printf("%s",m.c_str());
   fflush(stdout);
}

// ========================
// Write Application Header
// ========================
void WriteAppHeader(){
  WriteMessage("\n");
  WriteMessage("------------------------\n");
  WriteMessage("FEM Morphing Application\n");
  WriteMessage("Daniele E. Schiavazzi\n");
  WriteMessage("University of Notre Dame\n");
  WriteMessage("------------------------\n");
  WriteMessage("\n");
}

// ======================
// Write Application Help
// ======================
void WriteAppHelp(){
  WriteMessage("Usage: feMorph [OPTIONS] -f Input/Node File -o Output/Element File\n");
  WriteMessage("Where [OPTIONS] can be one of the following:\n");
  WriteMessage("-n               Normal Execution. Read an input file for morphing and creates\n");
  WriteMessage("                 the CVPre files for the SV presolver.\n");
  WriteMessage("-c               Translation to CVPre. Translates two files with node coordinates\n");
  WriteMessage("                 and element connectivity to CVPre files.\n");
  WriteMessage("-e               Mesh Quality Extraction. Creates quality statistics\n");
  WriteMessage("                 for tetrahedral meshes.\n");
  WriteMessage("-l               Face List Matching. Maps two list of VTK face triangulations.\n");
  WriteMessage("-s               Convert Surface triangulation into a mesh by using the tetgen mesher.\n");
  WriteMessage("-x               Compute model expectations from weighted model list.\n");
  WriteMessage("-t               Defile tolerance for face list matching function.\n");
  WriteMessage("-d               Activates debug mode.\n");
  WriteMessage("-y               Runs explicit VMS fluid solver on model.txt.\n");
  WriteMessage("-h               This help.\n");
}

// Extral product between vectors
void Do3DExternalProduct(double* v1,double* v2,double* v3){
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v2[2]*v1[0];
  v3[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

// Eval Internal Product
double Do3DInternalProduct(double* v1, double* v2){
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// Mixed Product
double Do3DMixedProduct(double* v1, double* v2, double* v3){
  double v12[3];
  Do3DExternalProduct(v1,v2,v12);
  return Do3DInternalProduct(v12,v3);
}

// =================================
// SORT INTEGER VECTOR - BUBBLE SORT
// =================================
void bubbleSortIntVector(std::vector<int>& nodes){
  int first;
  int second;
  int temp;
  for(size_t loopA=0;loopA<nodes.size();loopA++){
    for(size_t loopB=0;loopB<nodes.size();loopB++){
      first = nodes[loopA];
      second = nodes[loopB];
      if(first>second){
        temp = nodes[loopA];
        nodes[loopA] = nodes[loopB];
        nodes[loopB] = temp;
      }
    }
  }
}

// ===================
// Normalize 3D Vector
// ===================
void Normalize3DVector(double* vector){
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

// =======================
// Allocate and Initialize
// =======================
void matZeros(femDoubleMat& mat,int row_size,int column_size){
  mat.resize(row_size);
  for(int loopA=0;loopA<row_size;loopA++){
    mat[loopA].resize(column_size);
    for(int loopB=0;loopB<column_size;loopB++){
      mat[loopA][loopB] = 0.0;
    }
  }
}

// =================================================
// General Rotation Around a Vector with quaternions
// The Angle is in degrees
// =================================================
void Rotate3DVectorAroundAxis(double* vec, double angle,double* axis){
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
void ExportStenosisBoxToVTK(std::string fileName, std::vector<femNode*> steNodeList){
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
void ExportReferenceToVTKLegacy(femInputData* data, std::string fileName){
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
void GetAverageNodeFromList(std::vector<femNode*> &nodeList,double* stenosisBoxCenter){
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
int trunc(double value){
  int result = std::floor(std::fabs(value));
  return (value < 0.0) ? -result : result;
}

// =================================
// Do the Euclidean Norm of a Vector
// =================================
double DoEucNorm(double* vec){
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

// ============================
// COMPUTE TWO-NORM OF A VECTOR
// ============================
double DoEucNorm(const femDoubleVec& vec){
  double res = 0.0;
  for(ulint loopA=0;loopA<vec.size();loopA++){
    res += vec[loopA] * vec[loopA];
  }
  return sqrt(res);
}

// ====================================================
// Find intersection between coordinate plane and point
// ====================================================
void findPointsToCoordinate1PlaneIntersection(double planeCoord, double* coords1, double* coords2, std::vector<femNode*> &intNodeList){
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
void evalCoordinatePlaneIntersections(double currSliceCoord, double *coords1, double* coords2, double* coords3, std::vector<femNode*> &intNodeList){
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
void PlotSlicesToVTK(std::string fileName, std::vector<femModelSlice*> &slices){

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

// ====================================================================
// Make list compact by eliminating double entries: NOTE: order changed
// ====================================================================
void MakeCompactList(const std::vector<int>& list, std::vector<int> &compactList){
  
  // Check list size
  if(list.size() == 0){
    throw femException("Empty list in femUtils::MakeCompactList.\n");
  }

  // Sort list
  std::vector<int> inputList(list);
  std::sort(inputList.begin(),inputList.end());

  // Create compact copy
  int count = 0;
  int currentNode = 0;
  int nextNode = 0;
  while(count<(int)inputList.size()-1){
    // Get Nodes
    currentNode = inputList[count];
    nextNode = inputList[count+1];
    // If different then print
    if(currentNode != nextNode){
      // Print to file
      compactList.push_back(currentNode);
    }
    // Increment Counter
    count++;
  }
  // Print Last
  compactList.push_back(inputList[count]);
}

// ==============================================
// Find inverse mapping defined by integer vector
// ==============================================
int findVectorInverseMapping(int value, const std::vector<int>& list){
  int count = 0;
  bool found = false;
  while((!found)&&(count<(int)list.size())){
    found = (value == list[count]);
    // Update Count
    if(!found){
      count++;
    }
  }
  if(!found){
    throw femException("Internal: could not find mapping in findVectorInverseMapping.\n");
  }
  // Return
  return count;
}

// ==========================================
// Create Inverse Node Numbering Relationship
// ==========================================
void MakeInverseList(std::vector<int> list,int totalNodes,std::vector<int> &inverseList){
  // Allocate and initialize elements of the inverse list
  inverseList.resize(totalNodes);
  // Init
  for(int loopA=0;loopA<totalNodes;loopA++){
    inverseList[loopA] = -1;
  }
  // Assign Inverse Relationshop
  int currentidx = 0;
  for(unsigned int loopA=0;loopA<list.size();loopA++){
    currentidx = list[loopA];
    inverseList[currentidx] = loopA;
  }
}

// ==========================
// Write Vector Graph To File
// ==========================
void WriteGraphToFile(std::string fileName, int vecSize, std::vector<double> &vecX, std::vector<double> &vecY){
  // Open Output File
    FILE* outFile;
    outFile = fopen(fileName.c_str(),"w");
    // Write Header
  for(int loopA=0;loopA<vecSize;loopA++){
    fprintf(outFile,"%d %e %e\n",loopA,vecX[loopA],vecY[loopA]);
  }
    // Close Output file
    fclose(outFile);
}

// ====================
// WRITE VECTOR TO FILE
// ====================
void writeVectorToFile(std::string fileName, const femDoubleVec& vec){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<vec.size();loopA++){
    fprintf(outFile,"%25.10e\n",vec[loopA]);
  }
  // Close Output file
  fclose(outFile);
}

// ====================
// WRITE MATRIX TO FILE
// ====================
void writeMatrixToFile(std::string fileName, const femDoubleMat& mat){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  for(int loopA=0;loopA<mat.size();loopA++){
    for(int loopB=0;loopB<mat[loopA].size();loopB++){
      fprintf(outFile,"%5.3e ",mat[loopA][loopB]);
    }
    fprintf(outFile,"\n");
  }
  // Close Output file
  fclose(outFile);
}

// ======================
// PRINT MATRIX TO SCREEN
// ======================
void printMatrix(const femDoubleMat& mat){
  printf("\n");
  for(int loopA=0;loopA<mat.size();loopA++){
    for(int loopB=0;loopB<mat[loopA].size();loopB++){
      if(loopB == mat[loopA].size()-1){
        printf("%8.3e\n",mat[loopA][loopB]);          
      }else{
        printf("%8.3e , ",mat[loopA][loopB]);
      }      
    }
  }
}

// ======================
// PRINT VECTOR TO SCREEN
// ======================
void printVector(const femDoubleVec& vec){
  printf("\n");
  for(int loopA=0;loopA<vec.size();loopA++){
    printf("%8.3e\n",vec[loopA]);
  }
}

// =====================================
// CHECK IF A POINT IS INSIDE THE LIMITS
// =====================================
bool isInsideLimits(double* centroid,double* limitBox){
  bool isInsideX = (centroid[0]>=limitBox[0])&&(centroid[0]<=limitBox[1]);
  bool isInsideY = (centroid[1]>=limitBox[2])&&(centroid[1]<=limitBox[3]);
  bool isInsideZ = (centroid[2]>=limitBox[4])&&(centroid[2]<=limitBox[5]);
  if((isInsideX)&&(isInsideY)&&(isInsideZ)){
    return true;
  }else{
    return false;
  }
}

// =================================
// CHECK IF TWO BOXES ARE COMPATIBLE
// =================================
bool AreCompatibleBoxes(double* firstBox, double* secondBox,double tolerance){
  // Eval the Two diagonals
  double firstDiagonal = sqrt((firstBox[1]-firstBox[0])*(firstBox[1]-firstBox[0]) + (firstBox[3]-firstBox[2])*(firstBox[3]-firstBox[2]) + (firstBox[5]-firstBox[4])*(firstBox[5]-firstBox[4]));
  double secondDiagonal = sqrt((secondBox[1]-secondBox[0])*(secondBox[1]-secondBox[0]) + (secondBox[3]-secondBox[2])*(secondBox[3]-secondBox[2]) + (secondBox[5]-secondBox[4])*(secondBox[5]-secondBox[4]));
  // If diagonal are zero return
  if((fabs(firstDiagonal)<kMathZero)||(fabs(secondDiagonal)<kMathZero)){
    return false;
  }
  // If the diagonal is different than return
  if((fabs(firstDiagonal-secondDiagonal)/((double)firstDiagonal))>0.2){
    return false;
  }
  // Eval the 2-norm in the difference on box locations
  double boxDistance = sqrt((firstBox[0]-secondBox[0])*(firstBox[0]-secondBox[0]) +
                            (firstBox[1]-secondBox[1])*(firstBox[1]-secondBox[1]) +
                            (firstBox[2]-secondBox[2])*(firstBox[2]-secondBox[2]) +
                            (firstBox[3]-secondBox[3])*(firstBox[3]-secondBox[3]) +
                            (firstBox[4]-secondBox[4])*(firstBox[4]-secondBox[4]) +
                            (firstBox[5]-secondBox[5])*(firstBox[5]-secondBox[5]));
  // Make a choice based on distance and diagonal
  if(boxDistance<firstDiagonal*tolerance){
    return true;
  }else{
    return false;
  }
}

void ReadListFromFile(std::string fileName, std::vector<std::string> &fileList1){
  // Declare input File
  std::ifstream infile;
  infile.open(fileName.c_str());

  // Clear File List
  fileList1.clear();

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Add to list
    fileList1.push_back(buffer);
  }

  // Close File
  infile.close();
}

// ================================================
// READ WEIGHTED FILE LIST - FILE + WEIGHT (DOUBLE)
// ================================================
void ReadWeightedFileList(std::string fileName, std::vector<femWeightedFileName*> &fileList){
  // Declare
  double currWeight1 = 0.0;
  double currWeight2 = 0.0;

  // Declare input File
  std::ifstream infile;
  infile.open(fileName.c_str());

  // Tokenized Strings
  std::vector<std::string> tokenizedString;

  // Clear File List
  fileList.clear();

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Read File Name and Associated Weight
    if(tokenizedString.size() > 1){
      if(tokenizedString.size() == 2){
        currWeight1 = atof(tokenizedString[1].c_str());
        currWeight2 = currWeight1;
      }else if(tokenizedString.size() == 3){
        currWeight1 = atof(tokenizedString[1].c_str());
        currWeight2 = atof(tokenizedString[2].c_str());
      }else{
        throw femException("Invalid Weighted File Format.\n");
      }
      // Add to list
      femWeightedFileName* wFile = new femWeightedFileName(tokenizedString[0],currWeight1,currWeight2);
      fileList.push_back(wFile);
    }else{
      throw femException("Invalid Weighted File Format.\n");
    }
  }
  // Close File
  infile.close();
}

// ================================
// COMPUTE FROBENIUS NORM OF MATRIX
// ================================
double getMatrixNorm(const femDoubleMat& mat){
  double norm = 0.0;
  for(int loopA=0;loopA<mat.size();loopA++){
    for(int loopB=0;loopB<mat[loopA].size();loopB++){
      norm += mat[loopA][loopB] * mat[loopA][loopB];
    }
  }
  return sqrt(norm);
}


// Extract File Name From String
std::string extractFileName(std::string fullPath){
  std::vector<std::string> tokenizedString;
  // Trim String
  boost::trim(fullPath);
  // Tokenize String
  boost::split(tokenizedString,fullPath, boost::is_any_of("/"), boost::token_compress_on);
  // Return Last
  return tokenizedString[tokenizedString.size()-1];
}

// INVERT 3x3 MATRIX
void invert3x3Matrix(const femDoubleMat& m,femDoubleMat& invMat, double& detJ){

  double m00 = m[0][0];
  double m01 = m[0][1];
  double m02 = m[0][2];
  double m10 = m[1][0];
  double m11 = m[1][1];
  double m12 = m[1][2];
  double m20 = m[2][0];
  double m21 = m[2][1];
  double m22 = m[2][2];

  // computes the inverse of a matrix m
  detJ = m00 * (m11 * m22 - m21 * m12) -
         m01 * (m10 * m22 - m12 * m20) +
         m02 * (m10 * m21 - m11 * m20);

  double invdet = 1 / detJ;

  invMat.resize(kDims);
  for(int loopA=0;loopA<kDims;loopA++){
    invMat[loopA].resize(kDims);
  }

  // Matrix33d minv; // inverse of matrix m
  invMat[0][0] = (m11 * m22 - m21 * m12) * invdet;
  invMat[0][1] = (m02 * m21 - m01 * m22) * invdet;
  invMat[0][2] = (m01 * m12 - m02 * m11) * invdet;
  invMat[1][0] = (m12 * m20 - m10 * m22) * invdet;
  invMat[1][1] = (m00 * m22 - m02 * m20) * invdet;
  invMat[1][2] = (m10 * m02 - m00 * m12) * invdet;
  invMat[2][0] = (m10 * m21 - m20 * m11) * invdet;
  invMat[2][1] = (m20 * m01 - m00 * m21) * invdet;
  invMat[2][2] = (m00 * m11 - m10 * m01) * invdet;

  // double A11 = a22 * a33 - a23 * a32;
  // double A22 = a33 * a11 - a31 * a13;
  // double A33 = a11 * a22 - a12 * a21;
  // double A12 = a23 * a31 - a21 * a33;
  // double A23 = a31 * a12 - a32 * a11;
  // double A31 = a12 * a23 - a13 * a22;
  // double A21 = a32 * a13 - a12 * a33;
  // double A32 = a13 * a21 - a23 * a11;
  // double A13 = a21 * a22 - a31 * a22;
  // detJ = a11 * A11 + a12 * A21 + a13 * A31;
  // invMat[0][0] = A11/detJ;
  // invMat[0][1] = A12/detJ;
  // invMat[0][2] = A13/detJ;
  // invMat[1][0] = A21/detJ;
  // invMat[1][1] = A22/detJ;
  // invMat[1][2] = A23/detJ;
  // invMat[2][0] = A31/detJ;
  // invMat[2][1] = A32/detJ;
  // invMat[2][2] = A33/detJ;
}

// HYPERBOLIC COTANGENT
double coth(double value){
  double temp;
  temp = exp(value);
  return (temp + 1 / temp) / (temp - 1 / temp);
}

void invert3x3MatrixFor1DElements(const femDoubleMat& mat,femDoubleMat& invMat, double& detJ){
  invMat.resize(kDims);
  for(int loopA=0;loopA<kDims;loopA++){
    invMat[loopA].resize(kDims);
    for(int loopB=0;loopB<kDims;loopB++){
      invMat[loopA][loopB] = 0.0;
    }
  }
  detJ = 1.0;
  double derX = mat[0][0];
  double derY = mat[1][0];
  double derZ = mat[2][0];

  if(fabs(derX)>1.0e-8){
    invMat[0][0] = 1/derX;
    detJ *= derX;
  }
  if(fabs(derY)>1.0e-8){
    invMat[1][0] = 1/derY;
    detJ *= derY;
  }
  if(fabs(derZ)>1.0e-8){
    invMat[2][0] = 1/derZ;
    detJ *= derZ;
  }
}

string intToStr(int a){
  std::stringstream ss;
  ss << a;
  return ss.str();
}

string floatToStr(float a){
  std::stringstream ss;
  ss << a;
  return ss.str();
}

femDoubleVec interpVelocity(const double time, const double time1, const femDoubleVec& vel1, 
                                                      const double time2, const femDoubleVec& vel2, 
                                                      int dims){
  // Assume 3D velocities
  // Assume time1 < time < time2
  femDoubleVec res(dims,0.0);
  for(ulint loopA=0;loopA<dims;loopA++){
    res[loopA] = vel1[loopA] + ((vel2[loopA] - vel1[loopA])/(time2 - time1))*(time-time1);
  }
  return res;
}

femDoubleVec interpTableData(double time, const femDoubleMat& table, int time_idx){

  femDoubleVec vel1(time_idx,0.0);
  femDoubleVec vel2(time_idx,0.0);

  // If time is greater than table: return last value as a constant
  if(time >= table[table.size()-1][time_idx]){
    for(ulint loopA=0;loopA<time_idx;loopA++){
      vel1[loopA] = table[table.size()-1][loopA];
    }
    return vel1;
  }
  
  // Otherwise look for values less than and greater than current time
  bool found = false;
  ulint count = 1;
  while((!found)&&(count<table.size())){
    found = (time >= table[count-1][time_idx]) && (time < table[count][time_idx]);
    if(!found){
      count++;
    }
  }
  // If not found than terminate
  if(!found){
    printf("Time: %f\n",time);
    for(ulint loopA=0;loopA<table.size();loopA++){
      printf("VX: %f ",table[loopA][0]);
      printf("VY: %f ",table[loopA][1]);
      printf("VZ: %f ",table[loopA][2]);
      printf("V Time: %f\n",table[loopA][3]);
    }
    throw femException("Could not find time in femUtils::interpTableData.\n");
  }
  // Get the two velocities
  for(ulint loopA=0;loopA<time_idx;loopA++){
    vel1[loopA] = table[count-1][loopA];
    vel2[loopA] = table[count][loopA];
  }
  // Return interpolated velocity
  femDoubleVec res = interpVelocity(time, table[count-1][time_idx], vel1, table[count][time_idx], vel2);
  return res;
}

double getMaxModule(const femDoubleVec& vec){
  double res = 0.0;
  for(ulint i=0;i<vec.size();i++){
    if(fabs(vec[i] > res)){
      res = fabs(vec[i]);
    }
  }
  return res;
}

double getMaxModule(const femDoubleMat& mat, int dim1, int dim2){
  // Compute the norm over dimension 2
  femDoubleVec temp;
  double sum = 0.0;
  for(ulint a=0;a<dim1;a++){
    sum = 0.0;
    for(ulint i=0;i<dim2;i++){
      sum += mat[a][i]*mat[a][i];
    }
    temp.push_back(sqrt(sum));
  }
  // Take the maximum over dimension 1
  double res = 0.0;
  for(ulint a=0;a<dim1;a++){
    if(temp[a] > res){
      res = temp[a];
    }
  }
  return res;
}

bool file_exists(const string& fileName){
    std::ifstream infile(fileName.c_str());
    return infile.good();
}


}

#endif // FEMUTILS_CPP

