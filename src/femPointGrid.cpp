#include "femPointGrid.h"

// Create a regular grid
femPointGrid::femPointGrid(double* locMinPoint,
                           int* locTotPoints,
                           double* locGridAxis_S,
                           double* locGridAxis_T,
                           double* locGridAxis_U,
                           const femIntVec& dispNodes,
                           const femDoubleMat& dispVals){

  // Assign Min Points
  minPoint[0] = locMinPoint[0];
  minPoint[1] = locMinPoint[1];
  minPoint[2] = locMinPoint[2];

  // Assign Size
  size[0] = femUtils::DoEucNorm(locGridAxis_S);
  size[1] = femUtils::DoEucNorm(locGridAxis_T);
  size[2] = femUtils::DoEucNorm(locGridAxis_U);

  // Assign Centerpoint
  centerPoint[0] = minPoint[0] + 0.5 * size[0];
  centerPoint[1] = minPoint[1] + 0.5 * size[1];
  centerPoint[2] = minPoint[2] + 0.5 * size[2];

  // Assign Totals
  totPoints[0] = locTotPoints[0];
  totPoints[1] = locTotPoints[1];
  totPoints[2] = locTotPoints[2];

  int totalPoints = totPoints[0] * totPoints[1] * totPoints[2];

  // Assign Reference Grid
  gridAxis[0][0] = locGridAxis_S[0];
  gridAxis[1][0] = locGridAxis_S[1];
  gridAxis[2][0] = locGridAxis_S[2];
  gridAxis[0][1] = locGridAxis_T[0];
  gridAxis[1][1] = locGridAxis_T[1];
  gridAxis[2][1] = locGridAxis_T[2];
  gridAxis[0][2] = locGridAxis_U[0];
  gridAxis[1][2] = locGridAxis_U[1];
  gridAxis[2][2] = locGridAxis_U[2];

  coords.resize(totalPoints);
  disps.resize(totalPoints);
  for(int loopA=0;loopA<totalPoints;loopA++){
    coords[loopA].resize(3);
    coords[loopA][0] = 0.0;
    coords[loopA][1] = 0.0;
    coords[loopA][2] = 0.0;
    disps[loopA].resize(3);
    disps[loopA][0] = 0.0;
    disps[loopA][1] = 0.0;
    disps[loopA][2] = 0.0;
  }

  int count = 0;
  for(int loopA=0;loopA<totPoints[2];loopA++){
    for(int loopB=0;loopB<totPoints[1];loopB++){
      for(int loopC=0;loopC<totPoints[0];loopC++){
        coords[count][0] = minPoint[0];
        coords[count][1] = minPoint[1];
        coords[count][2] = minPoint[2];

        coords[count][0] += loopC/(double)(totPoints[0]-1) * gridAxis[0][0];
        coords[count][1] += loopC/(double)(totPoints[0]-1) * gridAxis[1][0];
        coords[count][2] += loopC/(double)(totPoints[0]-1) * gridAxis[2][0];

        coords[count][0] += loopB/(double)(totPoints[1]-1) * gridAxis[0][1];
        coords[count][1] += loopB/(double)(totPoints[1]-1) * gridAxis[1][1];
        coords[count][2] += loopB/(double)(totPoints[1]-1) * gridAxis[2][1];

        coords[count][0] += loopA/(double)(totPoints[2]-1) * gridAxis[0][2];
        coords[count][1] += loopA/(double)(totPoints[2]-1) * gridAxis[1][2];
        coords[count][2] += loopA/(double)(totPoints[2]-1) * gridAxis[2][2];

        count++;
      }
    }
  }

  double norm_S[3],norm_T[3],norm_U[3];
  norm_S[0] = gridAxis[0][0];
  norm_S[1] = gridAxis[1][0];
  norm_S[2] = gridAxis[2][0];
  norm_T[0] = gridAxis[0][1];
  norm_T[1] = gridAxis[1][1];
  norm_T[2] = gridAxis[2][1];
  norm_U[0] = gridAxis[0][2];
  norm_U[1] = gridAxis[1][2];
  norm_U[2] = gridAxis[2][2];
  femUtils::Normalize3DVector(norm_S);
  femUtils::Normalize3DVector(norm_T);
  femUtils::Normalize3DVector(norm_U);

  // Add Displacements along the Box axis
  for(int loopA=0;loopA<dispNodes.size();loopA++){

    disps[dispNodes[loopA]][0] += dispVals[loopA][0] * norm_S[0];
    disps[dispNodes[loopA]][1] += dispVals[loopA][0] * norm_S[1];
    disps[dispNodes[loopA]][2] += dispVals[loopA][0] * norm_S[2];

    disps[dispNodes[loopA]][0] += dispVals[loopA][1] * norm_T[0];
    disps[dispNodes[loopA]][1] += dispVals[loopA][1] * norm_T[1];
    disps[dispNodes[loopA]][2] += dispVals[loopA][1] * norm_T[2];

    disps[dispNodes[loopA]][0] += dispVals[loopA][2] * norm_U[0];
    disps[dispNodes[loopA]][1] += dispVals[loopA][2] * norm_U[1];
    disps[dispNodes[loopA]][2] += dispVals[loopA][2] * norm_U[2];

  }

}

// Empty Constructor
femPointGrid::~femPointGrid(){

}

// Simple Factorial
int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n-1)*n;
}

// Binomial Coefficient
double binomialCoeff(int n, int k){
  return factorial(n)/(factorial(n-k) * factorial(k));
}

// One dimensional Bernstein Basis
double eval1DBernsteinPolynomials(double normCoord,int i,int n){
  return binomialCoeff(n,i) * pow(normCoord,i) * pow(1.0 - normCoord,n-i);
}

void femPointGrid::evalTrilinearBernsteinPolynomials(double* currCoord,double* currResult){
  double currBBX,currBBY,currBBZ;
  int count = 0;
  // Loop on the coordinates
  for(int loopA=0;loopA<kDims;loopA++){
    currResult[loopA] = 0.0;
    // Loop on the three dimensions for the polynomial
    count = 0;
    for(int loopZ=0;loopZ<totPoints[2];loopZ++){
      for(int loopY=0;loopY<totPoints[1];loopY++){
        for(int loopX=0;loopX<totPoints[0];loopX++){
          // Eval the three Bernstein Basis
          currBBX = eval1DBernsteinPolynomials(currCoord[0],loopX,totPoints[0]-1);
          currBBY = eval1DBernsteinPolynomials(currCoord[1],loopY,totPoints[1]-1);
          currBBZ = eval1DBernsteinPolynomials(currCoord[2],loopZ,totPoints[2]-1);
          // Multivariate
          currResult[loopA] += currBBX * currBBY * currBBZ * disps[count][loopA];
          // Update Count
          count++;
        }
      }
    }
  }
}

// Check if the coords are valid
bool validBoxCoords(double s,double t,double u){
  bool valid_S = (s>=0.0)&&(s<=1.0);
  bool valid_T = (t>=0.0)&&(t<=1.0);
  bool valid_U = (u>=0.0)&&(u<=1.0);
  return valid_S && valid_T && valid_U;
}

// Evaluate the deformations at a set of point locations
void femPointGrid::morphPoints(const femDoubleMat& pointLocations, femDoubleMat& pointDeformations){
  double currCoord[3];
  double currResult[3];

  // Allocate Result
  pointDeformations.resize(pointLocations.size());
  for(int loopA=0;loopA<pointLocations.size();loopA++){
    pointDeformations[loopA].resize(3);
  }

  // Loop through the locations
  for(int loopA=0;loopA<pointLocations.size();loopA++){
    // Store the current coordinate
    currCoord[0] = pointLocations[loopA][0];
    currCoord[1] = pointLocations[loopA][1];
    currCoord[2] = pointLocations[loopA][2];

    // Check if the coordinates are inside the box
    if(validBoxCoords(currCoord[0],currCoord[1],currCoord[2])){
      // Evaluate the trivariate Bernstein Basis at this location
      evalTrilinearBernsteinPolynomials(currCoord,currResult);
      // Store the results
      pointDeformations[loopA][0] = currResult[0];
      pointDeformations[loopA][1] = currResult[1];
      pointDeformations[loopA][2] = currResult[2];
    }else{
      // No Change
      pointDeformations[loopA][0] = currCoord[0];
      pointDeformations[loopA][1] = currCoord[1];
      pointDeformations[loopA][2] = currCoord[2];
    }
    //printf("Disps: %f %f %f\n",currResult[0],currResult[1],currResult[2]);
    //getchar();
  }
}

// Get Box Coordinates
void femPointGrid::getBoxCoords(const femDoubleMat& spaceCoords, femDoubleMat& boxCoords){
  double S[3],T[3],U[3];
  double TU[3],SU[3],ST[3];
  double auxCoords[3];

  // Resize Result
  boxCoords.resize(spaceCoords.size());
  for(int loopA=0;loopA<spaceCoords.size();loopA++){
    boxCoords[loopA].resize(3);
  }

  // s = TxU.(X-X0)/TxU.S
  // t = SxU.(X-X0)/SxU.T
  // u = SxT.(X-X0)/SxT.U

  // S, T, U
  S[0] = gridAxis[0][0];
  S[1] = gridAxis[1][0];
  S[2] = gridAxis[2][0];
  T[0] = gridAxis[0][1];
  T[1] = gridAxis[1][1];
  T[2] = gridAxis[2][1];
  U[0] = gridAxis[0][2];
  U[1] = gridAxis[1][2];
  U[2] = gridAxis[2][2];

  // Evaluate the External Product
  femUtils::Do3DExternalProduct(T,U,TU);
  femUtils::Do3DExternalProduct(S,U,SU);
  femUtils::Do3DExternalProduct(S,T,ST);

  // Evaluate the Mixed Product
  double TUS = femUtils::Do3DMixedProduct(T,U,S);
  double SUT = femUtils::Do3DMixedProduct(S,U,T);
  double STU = femUtils::Do3DMixedProduct(S,T,U);

  for(int loopA=0;loopA<spaceCoords.size();loopA++){
    auxCoords[0] = spaceCoords[loopA][0] - minPoint[0];
    auxCoords[1] = spaceCoords[loopA][1] - minPoint[1];
    auxCoords[2] = spaceCoords[loopA][2] - minPoint[2];
    // Compute the inner product
    boxCoords[loopA][0] = femUtils::Do3DInternalProduct(TU,auxCoords)/TUS;
    boxCoords[loopA][1] = femUtils::Do3DInternalProduct(SU,auxCoords)/SUT;
    boxCoords[loopA][2] = femUtils::Do3DInternalProduct(ST,auxCoords)/STU;
  }
}

// Apply Free Form Deformation to a Model
void femPointGrid::morphModel(femModel* model){

  // Allocate and store nodes
  femDoubleMat meshNodesSpace;
  femDoubleMat meshNodesBox;
  femDoubleMat newMeshDisps;
  meshNodesSpace.resize(model->nodeList.size());
  for(int loopA=0;loopA<model->nodeList.size();loopA++){
    meshNodesSpace[loopA].resize(kDims);
    meshNodesSpace[loopA][0] = model->nodeList[loopA]->coords[0];
    meshNodesSpace[loopA][1] = model->nodeList[loopA]->coords[1];
    meshNodesSpace[loopA][2] = model->nodeList[loopA]->coords[2];
  }

  // Transform to box coords
  getBoxCoords(meshNodesSpace,meshNodesBox);

  // Morph all Nodes Using their Box Coordinates
  morphPoints(meshNodesBox, newMeshDisps);

  // Store Node Displacements
  for(int loopA=0;loopA<model->nodeList.size();loopA++){
    for(int loopB=0;loopB<kDims;loopB++){
      model->nodeList[loopA]->displacements[loopB] = newMeshDisps[loopA][loopB];
      model->nodeList[loopA]->displacements[loopB+kDims] = 0.0;
    }
  }
}

// Export Grid to VTK Legacy
void femPointGrid::ExportToVTKLegacy(std::string fileName){
  int totPts = totPoints[0] * totPoints[1] * totPoints[2];

  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Print Quantities
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Point Grid Print out\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");
  // List all the points with coordinates
  fprintf(outFile,"POINTS %d double\n",totPts);
  for(int loopA=0;loopA<totPts;loopA++){
    fprintf(outFile,"%e %e %e\n",coords[loopA][0],coords[loopA][1],coords[loopA][2]);
  }

  // Save Cells
  int node1,node2,node3,node4,node5,node6,node7,node8;
  int totCells = (totPoints[0]-1) * (totPoints[1]-1) * (totPoints[2]-1);
  fprintf(outFile,"CELLS %d %d\n",totCells,9*totCells);
  for(int loopA=0;loopA<totPoints[2]-1;loopA++){
    for(int loopB=0;loopB<totPoints[1]-1;loopB++){
      for(int loopC=0;loopC<totPoints[0]-1;loopC++){
        node1 = loopA * (totPoints[1] * totPoints[0]) + loopB * (totPoints[0]) + loopC;
        node2 = loopA * (totPoints[1] * totPoints[0]) + loopB * (totPoints[0]) + loopC + 1;
        node4 = loopA * (totPoints[1] * totPoints[0]) + (loopB + 1) * (totPoints[0]) + loopC;
        node3 = loopA * (totPoints[1] * totPoints[0]) + (loopB + 1) * (totPoints[0]) + loopC + 1;
        node5 = (loopA + 1) * (totPoints[1] * totPoints[0]) + loopB * (totPoints[0]) + loopC;
        node6 = (loopA + 1) * (totPoints[1] * totPoints[0]) + loopB * (totPoints[0]) + loopC + 1;
        node8 = (loopA + 1) * (totPoints[1] * totPoints[0]) + (loopB + 1) * (totPoints[0]) + loopC;
        node7 = (loopA + 1) * (totPoints[1] * totPoints[0]) + (loopB + 1) * (totPoints[0]) + loopC + 1;
        // Totals
        fprintf(outFile,"%d ",8);
        fprintf(outFile,"%d ",node1);
        fprintf(outFile,"%d ",node2);
        fprintf(outFile,"%d ",node3);
        fprintf(outFile,"%d ",node4);
        fprintf(outFile,"%d ",node5);
        fprintf(outFile,"%d ",node6);
        fprintf(outFile,"%d ",node7);
        fprintf(outFile,"%d\n",node8);
      }
    }
  }
  // Cell Types
  fprintf(outFile,"CELL_TYPES %d\n",totCells);
  for(int loopA=0;loopA<totCells;loopA++){
    fprintf(outFile,"%d\n",12);
  }
  // Control Points Displacements
  fprintf(outFile,"POINT_DATA %d\n",totPts);
  fprintf(outFile,"VECTORS DXYZ double\n");
  for(unsigned int loopA=0;loopA<totPts;loopA++){
    fprintf(outFile,"%e %e %e\n",disps[loopA][0],disps[loopA][1],disps[loopA][2]);
  }

  // Close File
  fclose(outFile);
  femUtils::WriteMessage(std::string("Point Grid Data Exported.\n"));
}


