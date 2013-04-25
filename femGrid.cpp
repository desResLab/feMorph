#include "femGrid.h"
#include "femModel.h"
#include "femConstants.h"

#include <math.h>

// Constructor for femGrid
femGrid::femGrid(femModel* model){
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
  double currCoord[3] = {0.0};
  int currNode = 0;
  int currIndex = 0;
  // TEMP: Grid with only corner nodes
  for(unsigned int loopA=0;loopA<model->elementList.size();loopA++){
    for(unsigned int loopB=0;loopB<model->elementList[loopA]->elementConnections.size();loopB++){
      //
      // Get Current Node
      currNode = model->elementList[loopA]->elementConnections[loopB];

      // Get Node Coord
      currCoord[0] = model->nodeList[currNode]->coords[0];
      currCoord[1] = model->nodeList[currNode]->coords[1];
      currCoord[2] = model->nodeList[currNode]->coords[2];

      // Get Index
      currIndex = ToIndexes(currCoord);

      if(currIndex > -1){
        bool found = false;
        unsigned int listCount = 0;
        while ((!found)&&(listCount<gridData[currIndex]->gridElementList.size())){
          // Check if found
          found = (gridData[currIndex]->gridElementList[listCount] == loopA);
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
    idx0 = (int)((coord[0] - gridLimits[0])/gridSpacing[0]);
  }
  // Idx1
  if (fabs(coord[1] - gridLimits[3])<kGridTol) {
    idx1 = (kGridSizeY-1);
  }else if (fabs(coord[1] - gridLimits[2])<kGridTol){
    idx1 = 0;
  }else{
    idx1 = (int)((coord[1] - gridLimits[2])/gridSpacing[1]);
  }
  // Idx2
  if (fabs(coord[2] - gridLimits[5])<kGridTol) {
    idx2 = (kGridSizeZ-1);
  }else if(fabs(coord[2] - gridLimits[4])<kGridTol){
    idx2 = 0;
  }else{
    idx2 = (int)((coord[2] - gridLimits[4])/gridSpacing[2]);
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

// Export Grid to VTK Legacy
void femGrid::ExportToVTKLegacy(std::string fileName){
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Print Quantities
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Grid Print out\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET STRUCTURED_POINTS\n");
  fprintf(outFile,"DIMENSIONS %d %d %d\n",kGridSizeX,kGridSizeY,kGridSizeZ);
  fprintf(outFile,"SPACING %e %e %e\n",gridSpacing[0],gridSpacing[1],gridSpacing[2]);
  fprintf(outFile,"ORIGIN %e %e %e\n",gridLimits[0],gridLimits[2],gridLimits[4]);
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
}

