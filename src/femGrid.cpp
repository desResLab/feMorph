#include "femGrid.h"
#include "femModel.h"
#include "femConstants.h"
#include "femUtils.h"

#include <math.h>

// Constructor for femGrid
femGrid::femGrid(femModel* model){
  // Write Message
  femUtils::WriteMessage(std::string("Creating Grid ..."));

  // Copy the modelBox from the Model
  for(int loopA=0;loopA<6;loopA++){
    gridLimits[loopA] = model->modelBox[loopA];
  }

  // Grid Spacing
  gridSpacing[0] = (gridLimits[1] - gridLimits[0])/(double(kGridSizeX));
  gridSpacing[1] = (gridLimits[3] - gridLimits[2])/(double(kGridSizeY));
  gridSpacing[2] = (gridLimits[5] - gridLimits[4])/(double(kGridSizeZ));

  // Allocate gridData
  int totalGridSize = kGridSizeX*kGridSizeY*kGridSizeZ;
  for(int loopA=0;loopA<totalGridSize;loopA++){
    femGridCell* gridCell = new femGridCell();
    gridData.push_back(gridCell);
  }

  // Fill Grid with element numbers
  double minCoord[3] = {0.0};
  double maxCoord[3] = {0.0};
  int minIndex[3] = {0};
  int maxIndex[3] = {0};
  int currIndex = 0;
  std::vector<femNode*> minMaxNodeList;
  // TEMP: Grid with only corner nodes
  for(unsigned int loopA=0;loopA<model->elementList.size();loopA++){

    // Create MinMax Node List
    minMaxNodeList.clear();
    model->elementList[loopA]->CreateMinMaxNodeList(model->nodeList,minMaxNodeList);

    // Get Min Node Coord
    minCoord[0] = minMaxNodeList[0]->coords[0];
    minCoord[1] = minMaxNodeList[0]->coords[1];
    minCoord[2] = minMaxNodeList[0]->coords[2];
    // Get Max Node Coord
    maxCoord[0] = minMaxNodeList[1]->coords[0];
    maxCoord[1] = minMaxNodeList[1]->coords[1];
    maxCoord[2] = minMaxNodeList[1]->coords[2];

    // Convert to maxMinIndex
    ToIndexArray(minCoord,minIndex);
    ToIndexArray(maxCoord,maxIndex);

    // Loop through the indexes
    for(int loopB=minIndex[0];loopB<(maxIndex[0]+1);loopB++){
      for(int loopC=minIndex[1];loopC<(maxIndex[1]+1);loopC++){
        for(int loopD=minIndex[2];loopD<(maxIndex[2]+1);loopD++){
          // Get Index Back
          currIndex = IndexArrayToIndex(loopB,loopC,loopD);
          //
          if(currIndex > -1){
            bool found = false;
            unsigned int listCount = 0;
            while ((!found)&&(listCount<gridData[currIndex]->gridElementList.size())){
              // Check if found
              found = (gridData[currIndex]->gridElementList[listCount] == (int)loopA);
              // Update ListCount
              if (!found){
                listCount++;
              }
            }
            // If not Found Insert
            if(!found){
              // Insert Element
              gridData[currIndex]->gridElementList.push_back(loopA);
            }
          }
        }
      }
    }
  }
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}

// Destructor
femGrid::~femGrid(){
}

// Convert Spatial Coordinates to Index
int femGrid::ToIndexes(double* coord){
  int idx0 = 0;
  int idx1 = 0;
  int idx2 = 0;
  // Idx0
  if (fabs(coord[0] - gridLimits[1])<kGridTol) {
    idx0 = (kGridSizeX-1);
  }else if (fabs(coord[0] - gridLimits[0])<kGridTol){
    idx0 = 0;
  }else{
    idx0 = femUtils::trunc((coord[0] - gridLimits[0])/gridSpacing[0]);
  }
  // Idx1
  if (fabs(coord[1] - gridLimits[3])<kGridTol) {
    idx1 = (kGridSizeY-1);
  }else if (fabs(coord[1] - gridLimits[2])<kGridTol){
    idx1 = 0;
  }else{
    idx1 = femUtils::trunc((coord[1] - gridLimits[2])/gridSpacing[1]);
  }
  // Idx2
  if (fabs(coord[2] - gridLimits[5])<kGridTol) {
    idx2 = (kGridSizeZ-1);
  }else if(fabs(coord[2] - gridLimits[4])<kGridTol){
    idx2 = 0;
  }else{
    idx2 = femUtils::trunc((coord[2] - gridLimits[4])/gridSpacing[2]);
  }
  // Check if Inside Limits
  if((idx0 > (kGridSizeX - 1))||(idx0 < 0)){
    return -1;
  }else if((idx1 > (kGridSizeY - 1))||(idx1 < 0)){
    return -1;
  }else if((idx2 > (kGridSizeZ - 1))||(idx2 < 0)){
    return -1;
  }else{
    return kGridSizeY*kGridSizeX*idx2 + kGridSizeX*idx1 + idx0;
  }
}

// =================================
// Convert Coordinate to Index Array
// =================================
void femGrid::ToIndexArray(double* coord,int* index){
  int idx0 = 0;
  int idx1 = 0;
  int idx2 = 0;
  // Idx0
  if (fabs(coord[0] - gridLimits[1])<kGridTol) {
    idx0 = (kGridSizeX-1);
  }else if (fabs(coord[0] - gridLimits[0])<kGridTol){
    idx0 = 0;
  }else{
    idx0 = femUtils::trunc((coord[0] - gridLimits[0])/gridSpacing[0]);
  }
  // Idx1
  if (fabs(coord[1] - gridLimits[3])<kGridTol) {
    idx1 = (kGridSizeY-1);
  }else if (fabs(coord[1] - gridLimits[2])<kGridTol){
    idx1 = 0;
  }else{
    idx1 = femUtils::trunc((coord[1] - gridLimits[2])/gridSpacing[1]);
  }
  // Idx2
  if (fabs(coord[2] - gridLimits[5])<kGridTol) {
    idx2 = (kGridSizeZ-1);
  }else if(fabs(coord[2] - gridLimits[4])<kGridTol){
    idx2 = 0;
  }else{
    idx2 = femUtils::trunc((coord[2] - gridLimits[4])/gridSpacing[2]);
  }
  // Save index Values
  index[0] = idx0;
  index[1] = idx1;
  index[2] = idx2;
}

// ================================
// From Index array to single index
// ================================
int femGrid::IndexArrayToIndex(int idx0,int idx1,int idx2){
  // Check if Inside Limits
  if((idx0 > (kGridSizeX - 1))||(idx0 < 0)){
    return -1;
  }else if((idx1 > (kGridSizeY - 1))||(idx1 < 0)){
    return -1;
  }else if((idx2 > (kGridSizeZ - 1))||(idx2 < 0)){
    return -1;
  }else{
    return kGridSizeY*kGridSizeX*idx2 + kGridSizeX*idx1 + idx0;
  }
}


// Export Grid to VTK Legacy
void femGrid::ExportToVTKLegacy(std::string fileName){
  // Write Message
  femUtils::WriteMessage(std::string("(debug) Exporting Grid to VTK ..."));

  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Print Quantities
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Grid Print out\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET STRUCTURED_POINTS\n");
  fprintf(outFile,"DIMENSIONS %d %d %d\n",kGridSizeX,kGridSizeY,kGridSizeZ);
  fprintf(outFile,"SPACING %e %e %e\n",gridSpacing[0],gridSpacing[1],gridSpacing[2]);
  fprintf(outFile,"ORIGIN %e %e %e\n",gridLimits[0]+0.5*gridSpacing[0],gridLimits[2]+0.5*gridSpacing[1],gridLimits[4]+0.5*gridSpacing[2]);
  fprintf(outFile,"POINT_DATA %d\n",kGridSizeX*kGridSizeY*kGridSizeZ);
  fprintf(outFile,"SCALARS Contain_Elements float 1\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");
  for(unsigned int loopA=0;loopA<kGridSizeX*kGridSizeY*kGridSizeZ;loopA++){
    if(gridData[loopA]->gridElementList.size()>0){
      fprintf(outFile,"%e\n",1.0);
    }else{
      fprintf(outFile,"%e\n",0.0);
    }
  }
  // Close File
  fclose(outFile);
  femUtils::WriteMessage(std::string("Done.\n"));
}

