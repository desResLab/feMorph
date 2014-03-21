#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <vector>
#include <stack>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>

#include "femModel.h"
#include "femGrid.h"
#include "femPoint.h"
#include "femInputData.h"
#include "femException.h"
#include "femConstants.h"
#include "femUtils.h"
#include "femStatistics.h"

using namespace std;

// Constructor
femModel::femModel()
{
  // Set modelBox
  // Min X
  modelBox[0] = std::numeric_limits<double>::max();
  // Max X
  modelBox[1] = -std::numeric_limits<double>::max();
  // Min Y
  modelBox[2] = std::numeric_limits<double>::max();
  // Max Y
  modelBox[3] = -std::numeric_limits<double>::max();
  // Min Z
  modelBox[4] = std::numeric_limits<double>::max();
  // Max Z
  modelBox[5] = -std::numeric_limits<double>::max();

  // Initialize Model Centre
  for(int loopA=0;loopA<3;loopA++){
    modelCentre[loopA] = 0.0;
  }
}

// Distructor
femModel::~femModel()
{
}

// =================
// Eval Model Centre
// =================
void femModel::EvalModelCentre(){
  modelCentre[0] = 0.5*(modelBox[0]+modelBox[1]);
  modelCentre[1] = 0.5*(modelBox[2]+modelBox[3]);
  modelCentre[2] = 0.5*(modelBox[4]+modelBox[5]);
}

// ======================================
// Get Node ID in List from its numbering
// ======================================
int femModel::GetNodeIDFromNumber(int number){
  bool found = false;
  uint nodeCount = 0;
  while ((!found)&&(nodeCount<nodeList.size())){
    // Check If Found
    found = (nodeList[nodeCount]->nodeNumber == number);

    // Update Node Counter
    if (!found){
      nodeCount++;
    }
  }
  if (!found){
      throw femException("Internal: Cannot Find Node Number.\n");
  } else return nodeCount;
}

// =================================================
// Read the whole model from LSDYNA-style input file
// =================================================
void femModel::ReadModelFromFile(std::string fileName){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    
    // Parse tokenized String
    if (boost::to_upper_copy(tokenizedString[1]) == std::string("NODE")){
      // Read Nodes

    }else if (boost::to_upper_copy(tokenizedString[1]) == std::string("TETRA10")){
      // Read Tetra 10 Element
    }else if (boost::to_upper_copy(tokenizedString[1]) == std::string("PROP")){
      // Read Element Property
    }
  }
  
  // Close File
  infile.close();
}

// ============
// Rotate Model
// ============
void femModel::RotateModel(double angle, double* axis){
  double dispVec[3] = {0.0};
  double rotVec[3] = {0.0};
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    // Rotate
    femUtils::Rotate3DVectorAroundAxis(nodeList[loopA]->coords,angle,axis);
    // Make local Copy
    for(int loopB=0;loopB<3;loopB++){
      dispVec[loopB] = nodeList[loopA]->displacements[loopB];
      rotVec[loopB] = nodeList[loopA]->displacements[loopB+3];
    }
    // Rotate Displacements and Rotations
    femUtils::Rotate3DVectorAroundAxis(dispVec,angle,axis);
    femUtils::Rotate3DVectorAroundAxis(rotVec,angle,axis);
    // Copy back
    for(int loopB=0;loopB<3;loopB++){
      nodeList[loopA]->displacements[loopB] = dispVec[loopB];
      nodeList[loopA]->displacements[loopB+3] = rotVec[loopB];
    }
  }
}

// ===================================
// Read the node coordinates from file
// ===================================
void femModel::ReadNodeCoordsFromFile(std::string fileName, bool skipFirstRow){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Write Message
  femUtils::WriteMessage(std::string("Reading Coordinates for file ")+fileName+std::string("..."));

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  // Initialize
  int currNodeNumber = 0;
  int lineCount = 0;
  double currCoords[3] = {0.0};
  double currDisps[6] = {0.0};
  // Skip First Row If Required
  if(skipFirstRow){
    std::getline(infile,buffer);
  }
  while (std::getline(infile,buffer)){
    // Increment line count
    lineCount++;
    // Trim String
    boost::trim(buffer);
    if((buffer != "")&&(buffer.at(0) != '#')){
      // Tokenize String
      boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

      // Read node Number
      currNodeNumber = atoi(tokenizedString[0].c_str());

      // Read Coordinates
      currCoords[0] = atof(tokenizedString[1].c_str());
      currCoords[1] = atof(tokenizedString[2].c_str());
      currCoords[2] = atof(tokenizedString[3].c_str());

      // Create a new Node
      femNode* newNode = new femNode(currNodeNumber,currCoords,currDisps);

      // Add to the node List
      nodeList.push_back(newNode);
    }
  }

  // Eval Model Box
  EvalModelBox();

  // Eval Model Centre
  EvalModelCentre();

  // Close File
  infile.close();

  femUtils::WriteMessage("Done.\n");
}

// =====================================
// Read Element Connectivities From File
// =====================================
void femModel::ReadElementConnectionsFromFile(std::string fileName, bool skipFirstRow){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Write Message
  femUtils::WriteMessage(std::string("Reading Connectivity for file ")+fileName+std::string("..."));

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  int totNodes = 0;
  int lineCount = 0;
  int currElementNumber = 0;
  // Skip First Row If Required
  if(skipFirstRow){
    std::getline(infile,buffer);
  }
  while (std::getline(infile,buffer)){
    // Increment Line Count
    lineCount++;
    // Trim String
    boost::trim(buffer);
    if((buffer != "")&&(buffer.at(0) != '#')){
      // Tokenize String
      boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

      // Read element Number: IMPORTANT they start from 1 in the file
      currElementNumber = atoi(tokenizedString[0].c_str()) - 1;

      // Read Element Connections
      totNodes = tokenizedString.size() - 1;
      int* nodeConnections = new int[totNodes];
      for (int loopA=0;loopA<totNodes;loopA++){
        // node numbering starts from 1
        nodeConnections[loopA] = atoi(tokenizedString[loopA+1].c_str()) - 1;
      }

      // Create a new Element: Fixed Property Number for the time being!!!
      femElement* newElement;
      if(totNodes == kTetra4Nodes){
        newElement = new femTetra4(currElementNumber,1,totNodes,nodeConnections);
      }else if(totNodes == kTetra10Nodes){
        newElement = new femTetra10(currElementNumber,1,totNodes,nodeConnections);
      }else{
        throw femException("Error: Invalid Element Type");
      }

      // Add to the node List
      elementList.push_back(newElement);
    }
  }

  // Test Print out
  //FILE* debugFile;
  //debugFile = fopen("testOut.txt","w");
  //for(int loopA=0;loopA<elementList.size();loopA++){
  //  fprintf(debugFile,"%d %d\n",loopA,elementList[loopA]->elementNumber);
  // }
  //fclose(debugFile);

  // Close File
  infile.close();

  femUtils::WriteMessage("Done.\n");
}

// =================================
// Read Node Displacements From File
// =================================
void femModel::ReadNodeDisplacementsFromFile(std::string fileName, bool readRotations){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Write Message
  femUtils::WriteMessage(std::string("Reading Displacements for file ")+fileName+std::string("..."));

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  int currNodeNumber = 0;
  int mapNodeNumber = 0;
  double currNodeDX = 0.0;
  double currNodeDY = 0.0;
  double currNodeDZ = 0.0;
  double currNodeRX = 0.0;
  double currNodeRY = 0.0;
  double currNodeRZ = 0.0;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

    // Read node Number
    currNodeNumber = atoi(tokenizedString[0].c_str());

    // Read Coordinates
    currNodeDX = atof(tokenizedString[1].c_str());
    currNodeDY = atof(tokenizedString[2].c_str());
    currNodeDZ = atof(tokenizedString[3].c_str());

    // Read rotations if required
    if (readRotations) {
      currNodeRX = atof(tokenizedString[4].c_str());
      currNodeRY = atof(tokenizedString[5].c_str());
      currNodeRZ = atof(tokenizedString[6].c_str());
    } else {
      currNodeRX = 0.0;
      currNodeRY = 0.0;
      currNodeRZ = 0.0;
    }

    // Find the Corresponding node
    mapNodeNumber = GetNodeIDFromNumber(currNodeNumber);

    // Add to the node List
    nodeList[mapNodeNumber]->setDisplacements(currNodeDX,currNodeDY,currNodeDZ,currNodeRX,currNodeRY,currNodeRZ);
  }

  // Close File
  infile.close();

  femUtils::WriteMessage("Done.\n");
}

// ====================
// Write Coords To File
// ====================
void femModel::WriteNodeCoordsToFile(double dispFactor, std::string fileName){
  // Write Message
  femUtils::WriteMessage(std::string("Writing Node Coords to File ")+fileName+std::string("..."));
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
      fprintf(outFile,"%d %16.8e %16.8e %16.8e\n",nodeList[loopA]->nodeNumber,
              nodeList[loopA]->coords[0] + nodeList[loopA]->displacements[0]*dispFactor,
              nodeList[loopA]->coords[1] + nodeList[loopA]->displacements[1]*dispFactor,
              nodeList[loopA]->coords[2] + nodeList[loopA]->displacements[2]*dispFactor);
  }
  // Close Output file
  fclose(outFile);
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}

// =================================
// Write Element Connections To File
// =================================
void femModel::WriteElementConnectionsToFile(std::string fileName){
  // Write Message
  femUtils::WriteMessage(std::string("Writing Element Connections to File ")+fileName+std::string("..."));
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Element Connectivities to File
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    fprintf(outFile,"%d ",elementList[loopA]->elementNumber + 1);
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementConnections.size();loopB++){
      fprintf(outFile,"%d ",elementList[loopA]->elementConnections[loopB] + 1);
    }
    fprintf(outFile,"\n");
  }
  // Close Output file
  fclose(outFile);
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}

// ===================================================
// Map Displacements between a main and mapping models
// ===================================================
void femModel::MapDisplacements(femModel* MappingModel,
                                femInputData* data,
                                double dispScaleFactor){
  // Initialize Mapping Coords
  double nodeCoords[3] = {0.0};
  double nodeDisps[3] = {0.0};  
  // Create New Grid
  femGrid* grid = new femGrid(MappingModel);

  // Check If Correct
  grid->ExportToVTKLegacy(std::string("grid.vtk"));

  // Write Message
  femUtils::WriteMessage(std::string("Mapping Displacements ..."));

  // Loop through the nodes in the main model
  int elementID = 0;
  int currIdx = 0;
  int progress = 0;
  int progressCounted = 0;
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    progress = (int)((loopA/double(nodeList.size()-1))*100.0);
    if (((progress % 10) == 0)&&((progress / 10) != progressCounted)){
      progressCounted = (progress / 10);
      femUtils::WriteMessage(std::string(boost::lexical_cast<std::string>(progress))+".");
    }
    //femUtils::WriteMessage(std::string(boost::lexical_cast<std::string>(loopA))+" ");
    // Store nodes
    nodeCoords[0] = nodeList[loopA]->coords[0];
    nodeCoords[1] = nodeList[loopA]->coords[1];
    nodeCoords[2] = nodeList[loopA]->coords[2];
    //MappingModel->elementList[1000]->evalElementCentroid(MappingModel->nodeList,nodeCoords);
    // Get Index
    currIdx = grid->ToIndexes(nodeCoords);
    // Find the enclosing element
    if(currIdx > -1){
      elementID = MappingModel->FindEnclosingElementWithGrid(0.0,nodeCoords,grid->gridData[currIdx]->gridElementList);
      //femUtils::WriteMessage("IN\n");
    }else{
      elementID = -1;
      //femUtils::WriteMessage("OUTSIDE\n");
    }
    if(elementID>-1){
      // Interpolate the displacements
      MappingModel->elementList[elementID]->InterpolateElementDisplacements(0.0,nodeCoords,MappingModel->nodeList,nodeDisps);
      // Save Displacements
      nodeList[loopA]->setDisplacements(nodeDisps[0],nodeDisps[1],nodeDisps[2],0.0,0.0,0.0);
    }else{
      // Save Zero Displacements
      nodeList[loopA]->setDisplacements(0.0,0.0,0.0,0.0,0.0,0.0);
    }
  }
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
  // Deallocate Grid
  delete grid;
}

// ==================================
// Find Face Given the Attached Nodes
// ==================================
int FindFace(std::vector<femFace*> &faceList,std::vector<int> &nodes){
  unsigned int count = 0;
  bool hasSameNodes = false;
  bool found = false;
  int innerCount = 0;
  while ((!found)&&(count<faceList.size())){
    hasSameNodes = true;
    innerCount = 0;
    while((hasSameNodes)&&(innerCount<(int)faceList[count]->faceNodes.size())){
      hasSameNodes = (hasSameNodes)&&(faceList[count]->faceNodes[innerCount] == nodes[innerCount]);
      // Update
      innerCount++;
    }
    // Assign Found
    found = hasSameNodes;
    // Update
    count++;
  }
  if(found){
    return faceList[count-1]->number;
  }else{
    return -1;
  }
}

// ==================================
// Find Edge Given the Attached Nodes
// ==================================
int FindEdge(std::vector<femEdge*> &edgeList,std::vector<int> &nodes){
  unsigned int count = 0;
  bool hasSameNodes = false;
  bool found = false;
  int innerCount = 0;
  while ((!found)&&(count<edgeList.size())){
    hasSameNodes = true;
    innerCount = 0;
    while((hasSameNodes)&&(innerCount<(int)edgeList[count]->edgeNodes.size())){
      hasSameNodes = (hasSameNodes)&&(edgeList[count]->edgeNodes[innerCount] == nodes[innerCount]);
      // Update
      innerCount++;
    }
    // Assign Found
    found = hasSameNodes;
    // Update
    count++;
  }
  if(found){
    return edgeList[count-1]->number;
  }else{
    return -1;
  }
}

// ======================
// Form Element-Face List
// ======================
// WARNING: WORKS ONLY FOR TETRA4 AND TETRA10!!!
void femModel::FormElementFaceList(){
  // Write Message
  femUtils::WriteMessage(std::string("Forming Face List..."));
  int totalFaces = 0;
  int progress = 0;
  int progressCounted = 0;
  std::vector<std::vector<femFace*>> firstNodefaceList;
  firstNodefaceList.resize((int)nodeList.size());
  std::vector<int> nodes;
  int faceID = -2;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    // Progress
    progress = (int)((loopA/double(elementList.size()-1))*100.0);
    if (((progress % 10) == 0)&&((progress / 10) != progressCounted)){
      progressCounted = (progress / 10);
      femUtils::WriteMessage(std::string(boost::lexical_cast<std::string>(progress))+".");
    }
    // Loop on the faces
    for(int loopB=0;loopB<kTetraFaces;loopB++){
      nodes.clear();
      nodes.push_back(elementList[loopA]->elementConnections[loopB % (kTetraFaces)]);
      nodes.push_back(elementList[loopA]->elementConnections[(loopB + 1) % (kTetraFaces)]);
      nodes.push_back(elementList[loopA]->elementConnections[(loopB + 2) % (kTetraFaces)]);
      // Sort Nodes
      std::sort(std::begin(nodes), std::end(nodes));

      // Check if the face is already there
      faceID = FindFace(firstNodefaceList[nodes[0]],nodes);
      if (faceID>-1){
        elementList[loopA]->elementFaces.push_back(faceID);
        faceList[faceID]->faceElements.push_back(loopA);
      }else{
        totalFaces++;
        femFace* face = new femFace();
        faceList.push_back(face);
        elementList[loopA]->elementFaces.push_back(totalFaces-1);
        faceList[totalFaces-1]->faceElements.push_back(loopA);
        for(unsigned int loopC=0;loopC<nodes.size();loopC++){
          faceList[totalFaces-1]->faceNodes.push_back(nodes[loopC]);
        }
        // Add to the list associated to the first node
        firstNodefaceList[nodes[0]].push_back(new femFace(totalFaces-1,nodes));
      }
    }
  }
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}

// Sorting Function
bool orderTuple (std::tuple<int,int,int> i,std::tuple<int,int,int> j){
  return (std::get<0>(i)<std::get<0>(j))&&(std::get<1>(i)<std::get<1>(j))&&(std::get<2>(i)<std::get<2>(j));
}

// ========================================
// Get neighbor element of minimum distance
// ========================================
void femModel::getNextElement(int currElement, double* nodeCoords, int &nextElement, double &nextDistance){
  std::vector<int> elFaces;
  std::vector<int> searchedElements;
  // Get the list of element faces
  for(unsigned int loopA=0;loopA<elementList[currElement]->elementFaces.size();loopA++){
    elFaces.push_back(elementList[currElement]->elementFaces[loopA]);
  }
  // Get the list of connected elements
  int otherElement = 0;
  for(unsigned int loopA=0;loopA<elFaces.size();loopA++){
    for(unsigned int loopB=0;loopB<faceList[elFaces[loopA]]->faceElements.size();loopB++){
      otherElement = faceList[elFaces[loopA]]->faceElements[loopB];
      if(otherElement!=currElement){
        searchedElements.push_back(otherElement);
      }
    }
  }
  // Find the element Of Minimum Distance
  nextDistance = std::numeric_limits<double>::max();
  double currDist = 0.0;
  for(unsigned int loopA=0;loopA<searchedElements.size();loopA++){
    currDist = elementList[searchedElements[loopA]]->evalPointToElementDistance(nodeCoords,nodeList);
    if(currDist<nextDistance){
      nextDistance = currDist;
      nextElement = searchedElements[loopA];
    }
  }
}

// =====================================
// Find the Enclosing element: Presearch
// =====================================
int femModel::FindEnclosingElementPre(double* nodeCoords){
  // Initialize
  int countMC = 0;
  int currElement = 0;
  double currDistance = 0.0;
  int nextElement = 0;
  double nextDistance = 0.0;
  bool endSearch = false;

  // Init min distance and min element
  double minDistance = std::numeric_limits<double>::max();
  int minElement = 0;

  // Loop on a Monte Carlo possible choice of nodes
  while (countMC<kMaxEnclosingMC){
    // Start from a random Element
    currElement = femUtils::GenerateUniformIntegers(0,elementList.size()-1);
    currDistance = elementList[currElement]->evalPointToElementDistance(nodeCoords,nodeList);
    // Search for a local miminum in the distance function
    endSearch = false;
    while(!endSearch){
      // Get Next Element
      getNextElement(currElement,nodeCoords,nextElement,nextDistance);
      endSearch = (nextDistance > currDistance);
      // Update
      if(!endSearch){
        currDistance = nextDistance;
        currElement = nextElement;
      }
    }
    // Store Minimum distance element
    if(currDistance<minDistance){
      minDistance = currDistance;
      minElement = currElement;
    }
    // Update
    countMC++;
  }
  // Return from function
  return minElement;
}

// =================================
// Find the Enclosing element: Check
// =================================
int femModel::FindEnclosingElementPost(double* nodeCoords, int startingElement){
  int currElement = startingElement;
  int previous = currElement;
  // Initialize is visited Property
  bool isVisited[elementList.size()];
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    isVisited[loopA] = false;
  }
  bool isVisited1 = false;
  bool isVisited2 = false;
  bool isVisited3 = false;
  bool isVisited4 = false;
  bool finished = false;
  int localFaceID = 0;
  bool isOnFace = false;
  bool isOnOppositeSide = false;
  int neighbor = 0;
  while(!finished){
    localFaceID = femUtils::GenerateUniformIntegers(0,kTetraFaces - 1);
    // Check two conditions
    elementList[currElement]->CheckRandomWalkingCriterion(localFaceID,nodeCoords,faceList,nodeList,isOnFace,isOnOppositeSide);
    // Check neighbor
    neighbor = elementList[currElement]->getAdjacentElement(localFaceID,faceList);
    // Check If visited
    if(neighbor > -1){
      isVisited1 = isVisited[neighbor];
    }
    if((!isOnFace)&&(isOnOppositeSide)&&(neighbor > -1)&&(!isVisited1)){
      previous = currElement;
      currElement = neighbor;
      isVisited[currElement] = true;
    }else{
      // Get next local face ID
      localFaceID = ((localFaceID + 1) % (kTetraFaces));
      // Check Conditions
      elementList[currElement]->CheckRandomWalkingCriterion(localFaceID,nodeCoords,faceList,nodeList,isOnFace,isOnOppositeSide);
      // Check neighbor
      neighbor = elementList[currElement]->getAdjacentElement(localFaceID,faceList);
      // Check If visited
      if(neighbor > -1){
        isVisited2 = isVisited[neighbor];
      }
      if((!isOnFace)&&(isOnOppositeSide)&&(neighbor > -1)&&(!isVisited2)){
        previous = currElement;
        currElement = neighbor;
        isVisited[currElement] = true;
      }else{
        // Get next local face ID
        localFaceID = ((localFaceID + 2) % (kTetraFaces));
        // Check Conditions
        elementList[currElement]->CheckRandomWalkingCriterion(localFaceID,nodeCoords,faceList,nodeList,isOnFace,isOnOppositeSide);
        // Check neighbor
        neighbor = elementList[currElement]->getAdjacentElement(localFaceID,faceList);
        // Check If visited
        if(neighbor > -1){
          isVisited3 = isVisited[neighbor];
        }
        if((!isOnFace)&&(isOnOppositeSide)&&(neighbor > -1)&&(!isVisited3)){
          previous = currElement;
          currElement = neighbor;
          isVisited[currElement] = true;
        }else{
          // Get next local face ID
          localFaceID = ((localFaceID + 3) % (kTetraFaces));
          // Check Conditions
          elementList[currElement]->CheckRandomWalkingCriterion(localFaceID,nodeCoords,faceList,nodeList,isOnFace,isOnOppositeSide);
          // Check neighbor
          neighbor = elementList[currElement]->getAdjacentElement(localFaceID,faceList);
          // Check If visited
          if(neighbor > -1){
            isVisited4 = isVisited[neighbor];
          }
          if((!isOnFace)&&(isOnOppositeSide)&&(neighbor > -1)&&(!isVisited4)){
            previous = currElement;
            currElement = neighbor;
            isVisited[currElement] = true;
          }else{
            finished = true;
          }
        }
      }
    }
  }
  // If has visited all
  if(isVisited1 && isVisited2 && isVisited3 && isVisited3){
    return -1;
  }else{
    return previous;
  }
}


// ==============
// Eval Model Box
// ==============
void femModel::EvalModelBox(){
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    // Min X
    if (nodeList[loopA]->coords[0]<modelBox[0]){
      modelBox[0] = nodeList[loopA]->coords[0];
    }
    // Max X
    if (nodeList[loopA]->coords[0]>modelBox[1]){
      modelBox[1] = nodeList[loopA]->coords[0];
    }
    // Min Y
    if (nodeList[loopA]->coords[1]<modelBox[2]){
      modelBox[2] = nodeList[loopA]->coords[1];
    }
    // Max Y
    if (nodeList[loopA]->coords[1]>modelBox[3]){
      modelBox[3] = nodeList[loopA]->coords[1];
    }
    // Min Z
    if (nodeList[loopA]->coords[2]<modelBox[4]){
      modelBox[4] = nodeList[loopA]->coords[2];
    }
    // Max Z
    if (nodeList[loopA]->coords[2]>modelBox[5]){
      modelBox[5] = nodeList[loopA]->coords[2];
    }
  }
}

// ===========================================
// Check if a node is outside the model limits
// ===========================================
bool femModel::IsInsideLimits(double* nodeCoords){
  bool isInside = false;
  isInside = ((nodeCoords[0]>=modelBox[0])&&(nodeCoords[0]<=modelBox[1])&&
             (nodeCoords[1]>=modelBox[2])&&(nodeCoords[1]<=modelBox[3])&&
             (nodeCoords[2]>=modelBox[4])&&(nodeCoords[2]<=modelBox[5]));
  // Return
  return isInside;
}

// ================
// Get Stenosis Box
// ================
void femModel::GetStenosisBox(femInputData* data, double* limRect){
  // Allocate Coords
  double newCoords[3] = {0.0};

  // Intialize Limits
  limRect[0] = -0.5 * data->stenosisLength;
  limRect[1] =  0.5 * data->stenosisLength;
  limRect[2] = std::numeric_limits<double>::max();
  limRect[3] = -std::numeric_limits<double>::max();
  limRect[4] = std::numeric_limits<double>::max();
  limRect[5] = -std::numeric_limits<double>::max();

  // Loop through all the nodes
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){

    // Transform Node Coords
    nodeList[loopA]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,newCoords,kDirect,kUndeformed,0.0);

    // If the first coord is within box than store the other two
    if ((newCoords[0]<(0.5*data->stenosisLength))&&(newCoords[0]>(-0.5*data->stenosisLength))&&
        (newCoords[1]<(kStenosisBoxFactor*data->stenosisLength))&&(newCoords[1]>(-kStenosisBoxFactor*data->stenosisLength))&&
        (newCoords[2]<(kStenosisBoxFactor*data->stenosisLength))&&(newCoords[2]>(-kStenosisBoxFactor*data->stenosisLength))){
      // Min 1
      if (newCoords[1]<limRect[2]){
         limRect[2] = newCoords[1];
      }
      // Max 1
      if (newCoords[1]>limRect[3]){
         limRect[3] = newCoords[1];
      }
      // Min 2
      if (newCoords[2]<limRect[4]){
         limRect[4] = newCoords[2];
      }
      // Max 2
      if (newCoords[2]>limRect[5]){
         limRect[5] = newCoords[2];
      }
    }
  }
}

// ============
// Store Limits
// ============
void updateLimits(double* newCoords, double* limRect){
  // Min 1
  if (newCoords[0]<limRect[0]){
    limRect[0] = newCoords[0];
  }
  // Min 1
  if (newCoords[0]>limRect[1]){
    limRect[1] = newCoords[0];
  }
  // Min 2
  if (newCoords[1]<limRect[2]){
    limRect[2] = newCoords[1];
  }
  // Max 2
  if (newCoords[1]>limRect[3]){
    limRect[3] = newCoords[1];
  }
  // Min 3
  if (newCoords[2]<limRect[4]){
    limRect[4] = newCoords[2];
  }
  // Max 3
  if (newCoords[2]>limRect[5]){
    limRect[5] = newCoords[2];
  }
}

// ==============================================================
// Find rotation axis and angles to match main and mapping models
// Note that the mapping model is in the standard reference frame
// ==============================================================
void getAxisAndRotations(double** refSystem, double &angle1, double &angle2, double* axis1, double* axis2){
  double mapVec1[3] = {0.0};
  double mapVec2[3] = {0.0};
  double secondVec[3] = {0.0};
  mapVec1[0] = 1.0;
  mapVec2[1] = 1.0;

  // Get first vector
  for(int loopA=0;loopA<3;loopA++){
    axis2[loopA] = refSystem[loopA][0];
    secondVec[loopA] = refSystem[loopA][1];
  }
  femUtils::Normalize3DVector(axis2);
  angle1 = acos(axis2[0])*(180.0/kPI);
  femUtils::Do3DExternalProduct(mapVec1,axis2,axis1);

  // Get Second Axis
  femUtils::Rotate3DVectorAroundAxis(mapVec2,angle1,axis1);
  femUtils::Normalize3DVector(mapVec2);
  femUtils::Normalize3DVector(secondVec);
  angle2 = acos(femUtils::Do3DInternalProduct(mapVec2,secondVec))*(180.0/kPI);
  femUtils::Do3DExternalProduct(mapVec2,secondVec,axis2);
}

// =========================================
// Transform node position and displacements
// =========================================
femModel* femModel::TransformModel(femInputData* data, double* stenosisBoxCenter, double* stenosisBox){

  // Create the new Model
  femModel* target = new femModel();

  // Copy everything but the nodes from the Current Model
  CopyElementsTo(target);
  CopyPropertyTo(target);

  // Intialize Limits
  double limRect[6];
  limRect[0] = std::numeric_limits<double>::max();
  limRect[1] = -std::numeric_limits<double>::max();
  limRect[2] = std::numeric_limits<double>::max();
  limRect[3] = -std::numeric_limits<double>::max();
  limRect[4] = std::numeric_limits<double>::max();
  limRect[5] = -std::numeric_limits<double>::max();

  // Tranlate Model with Centre in Origin
  double newCoords[3];
  double newDisps[6];
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){

    // Transform Node Coords
    for(int loopB=0;loopB<3;loopB++){
      newCoords[loopB] = nodeList[loopA]->coords[loopB] - modelCentre[loopB];
    }

    // Transform Displacements: They do not need to be rotated!
    for(int loopB=0;loopB<6;loopB++){
      newDisps[loopB] = nodeList[loopA]->displacements[loopB];
    }

    // Add to the Target Node List
    femNode* node = new femNode(loopA,newCoords,newDisps);
    target->nodeList.push_back(node);

    // Store Limits
    updateLimits(newCoords,limRect);
  }

  // Find Scale Factors
  double scaleFactor[3];
  for(int loopA=0;loopA<3;loopA++){
    scaleFactor[loopA] = ((stenosisBox[loopA*2]-stenosisBox[loopA*2+1])/double(limRect[loopA*2]-limRect[loopA*2+1]));
  }
  // Make same transverse scaling factors
  double maxTransverseScaling = max(scaleFactor[1],scaleFactor[2]);
  scaleFactor[1] = kTransverseScalingFactor*maxTransverseScaling;
  scaleFactor[2] = kTransverseScalingFactor*maxTransverseScaling;

  // Scale the model around the origin
  for(unsigned int loopA=0;loopA<target->nodeList.size();loopA++){
    for(int loopB=0;loopB<3;loopB++){
      // Scale and translate it to the new location
      target->nodeList[loopA]->coords[loopB] = target->nodeList[loopA]->coords[loopB]*scaleFactor[loopB];
    }
  }

  // Find rotation axis and angles
  double angle1 = 0.0;
  double angle2 = 0.0;
  double axis1[3] = {0.0};
  double axis2[3] = {0.0};
  getAxisAndRotations(data->mainModelRefSystem,angle1,angle2,axis1,axis2);

  // Rotate Model: First Rotation
  target->RotateModel(angle1,axis1);

  // Rotate Model: Second Rotation
  target->RotateModel(angle2,axis2);

  // Add Origin
  for(unsigned int loopA=0;loopA<target->nodeList.size();loopA++){
    for(int loopB=0;loopB<3;loopB++){
      target->nodeList[loopA]->coords[loopB] += stenosisBoxCenter[loopB];
    }
  }

  // Get Box and Centre
  target->EvalModelBox();
  target->EvalModelCentre();

  // Create Face List
  target->FormElementFaceList();

  // Return Transformed Model
  return target;
}

// ============================
// Copy elements to other model
// ============================
void femModel::CopyElementsTo(femModel* otherModel){
  femElement* newElement;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    if(elementList[loopA]->elementConnections.size() == 4){
      newElement = new femTetra4(elementList[loopA]);
    }else if(elementList[loopA]->elementConnections.size() == 10){
      newElement = new femTetra10(elementList[loopA]);
    }else{
      throw femException("Internal: Element not supported.");
    }
    // Put in list
    otherModel->elementList.push_back(newElement);
  }
}

// =========================
// Copy faces to other model
// =========================
void femModel::CopyFacesTo(femModel* otherModel){
  for(unsigned int loopA=0;loopA<faceList.size();loopA++){
    femFace* newFace = new femFace(faceList[loopA]);
    otherModel->faceList.push_back(newFace);
  }
}

// ============================
// Copy Property to other model
// ============================
void femModel::CopyPropertyTo(femModel* otherModel){
  for(unsigned int loopA=0;loopA<propList.size();loopA++){
    femProperty* newProp = new femProperty(propList[loopA]);
    otherModel->propList.push_back(newProp);
  }
}

// ==========================
// Export Model to VTK Legacy
// ==========================
void femModel::ExportToVTKLegacy(std::string fileName){
  // Write Message
  femUtils::WriteMessage(std::string("(Debug) Exporting Model to VTK file ")+fileName+std::string("..."));
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Wrtie Header
  fprintf(outFile,"# vtk DataFile Version 2.0\n");
  fprintf(outFile,"Model Exported from feMorph\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET UNSTRUCTURED_GRID\n");
  // Write Points
  fprintf(outFile,"POINTS %d float\n",int(nodeList.size()));
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    fprintf(outFile,"%e %e %e\n",nodeList[loopA]->coords[0],nodeList[loopA]->coords[1],nodeList[loopA]->coords[2]);
  }
  // Count the size of the cell list
  int cellListSize = 0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    cellListSize += elementList[loopA]->elementConnections.size() + 1;
  }
  // Write CELLS header
  fprintf(outFile,"CELLS %d %d\n",int(elementList.size()),cellListSize);
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    fprintf(outFile,"%d ",int(elementList[loopA]->elementConnections.size()));
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementConnections.size();loopB++){
      fprintf(outFile,"%d ",int(elementList[loopA]->elementConnections[loopB]));
    }
    fprintf(outFile,"\n");
  }
  // Write Cells Type Header
  fprintf(outFile,"CELL_TYPES %d\n",int(elementList.size()));
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    if(elementList[loopA]->elementConnections.size() == 4){
      fprintf(outFile,"%d\n",10);
    }else if(elementList[loopA]->elementConnections.size() == 10){
      fprintf(outFile,"%d\n",24);
    }else if(elementList[loopA]->elementConnections.size() == 3){
      fprintf(outFile,"%d\n",5);
    }else{
      fclose(outFile);
      throw femException("Error: Invalid element to Export.");
    }
  }
  // Point Data Header
  fprintf(outFile,"POINT_DATA %d\n",(int)nodeList.size());
  // Save Displacements DX,DY,DZ as vectors
  fprintf(outFile,"VECTORS DXYZ float\n");
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    fprintf(outFile,"%e %e %e\n",nodeList[loopA]->displacements[0],nodeList[loopA]->displacements[1],nodeList[loopA]->displacements[2]);
  }
  // Save properties as scalars
  fprintf(outFile,"CELL_DATA %d\n",int(elementList.size()));
  fprintf(outFile,"SCALARS cell_scalars int 1\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    fprintf(outFile,"%d\n",(int)elementList[loopA]->propertyNumber);
  }
  // Close Output file
  fclose(outFile);
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}


// =============================
// Create Node List for Stenosis
// =============================
void femModel::createStenosisNodeList(double* stenosisBox, femInputData* data, std::vector<femNode*> &steNodeList){
  double currDisps[3] = {0.0};
  double coord[3] = {0.0};
  double boxCoords[3] = {0.0};
  // Axes
  double axis1[3] = {0.0};
  double axis2[3] = {0.0};
  double axis3[3] = {0.0};
  // Axis 1
  axis1[0] = data->mainModelRefSystem[0][0];
  axis1[1] = data->mainModelRefSystem[1][0];
  axis1[2] = data->mainModelRefSystem[2][0];
  femUtils::Normalize3DVector(axis1);
  // Axis 2
  axis2[0] = data->mainModelRefSystem[0][1];
  axis2[1] = data->mainModelRefSystem[1][1];
  axis2[2] = data->mainModelRefSystem[2][1];
  femUtils::Normalize3DVector(axis2);
  // Axis 3
  axis3[0] = data->mainModelRefSystem[0][2];
  axis3[1] = data->mainModelRefSystem[1][2];
  axis3[2] = data->mainModelRefSystem[2][2];
  femUtils::Normalize3DVector(axis3);
  // Normalize Axis
  // Write Coords for every point
  for(int loopA=0;loopA<8;loopA++){
    switch(loopA){
      case 0:
        coord[0] = stenosisBox[0];
        coord[1] = stenosisBox[2];
        coord[2] = stenosisBox[4];
        break;
      case 1:
        coord[0] = stenosisBox[1];
        coord[1] = stenosisBox[2];
        coord[2] = stenosisBox[4];
        break;
      case 2:
        coord[0] = stenosisBox[0];
        coord[1] = stenosisBox[3];
        coord[2] = stenosisBox[4];
        break;
      case 3:
        coord[0] = stenosisBox[1];
        coord[1] = stenosisBox[3];
        coord[2] = stenosisBox[4];
        break;
      case 4:
        coord[0] = stenosisBox[0];
        coord[1] = stenosisBox[2];
        coord[2] = stenosisBox[5];
        break;
      case 5:
        coord[0] = stenosisBox[1];
        coord[1] = stenosisBox[2];
        coord[2] = stenosisBox[5];
        break;
      case 6:
        coord[0] = stenosisBox[0];
        coord[1] = stenosisBox[3];
        coord[2] = stenosisBox[5];
        break;
      case 7:
        coord[0] = stenosisBox[1];
        coord[1] = stenosisBox[3];
        coord[2] = stenosisBox[5];
        break;
    }
    // Add the coord and Points
    boxCoords[0] = data->mainModelOrigin[0] + coord[0]*axis1[0] + coord[1]*axis2[0] + coord[2]*axis3[0];
    boxCoords[1] = data->mainModelOrigin[1] + coord[0]*axis1[1] + coord[1]*axis2[1] + coord[2]*axis3[1];
    boxCoords[2] = data->mainModelOrigin[2] + coord[0]*axis1[2] + coord[1]*axis2[2] + coord[2]*axis3[2];
    // Add to Node List
    femNode* newNode = new femNode(loopA,boxCoords,currDisps);
    steNodeList.push_back(newNode);
  }
}

// =====================================
// Find Enclosing Element with Adjacency
// =====================================
int femModel::FindEnclosingElementWithAdj(double* nodeCoords){
  // Find the Enclosing element for the node
  bool isInside = IsInsideLimits(nodeCoords);
  if(!isInside){
    return -1;
  }else{
    // Find Element Containing a Given Node: Presearch
    int startingElement = FindEnclosingElementPre(nodeCoords);
    // Find Element Containing a Given Node: Walking Algortihm
    return FindEnclosingElementPost(nodeCoords,startingElement);
  }
}

// ================================
// Find Enclosing Element with Grid
// ================================
int femModel::FindEnclosingElementWithGrid(double dispFactor, double* nodeCoords, std::vector<int> &gridElementList){
  unsigned int count = 0;
  bool found = false;
  int currElement = 0;
  while ((!found)&&(count<gridElementList.size())){
    // Store element
    currElement = gridElementList[count];
    found = elementList[currElement]->isNodeInsideElement(dispFactor,nodeCoords,nodeList);
    // Update
    if(!found){
      count++;
    }
  }
  if(!found){
    return -1;
  }else{
    return currElement;
  }
}

// ================================
// Seek Stenosis levels iteratively
// ================================
double femModel::seekStenoticDisplacementFactor(femInputData* data, double targetStenosisLevel, bool debugMode){
  // Print Iteration Header
  printf("\n");
  printf("Picard Iterations for stenosis Level: %e\n",targetStenosisLevel);
  printf("%5s %20s %20s %20s\n","IT","DISP FACTOR","STENOSIS LEVEL","RESIDUAL");
  fflush(stdout);

  // Open the File if in debug mode
  FILE* debugFile;
  FILE* defFile;
  if (debugMode){
    debugFile = fopen("stenosisAreas.dat","a");
    defFile = fopen("stenosisDefs.dat","a");
  }
  // Print Stenosis Definitions header
  fprintf(defFile,"%15s %15s %15s %15s %15s %15s\n","Ref Diam","Ref Area","Old Diam","New Diam","Old Area","New Area");

  std::vector<femModelSlice*> slices;
  std::vector<double> sliceAreas;
  double stenosisDefs[6] = {0.0};
  bool converged = false;
  int currIt = 0;
  double firstDispFactor = 1.0;
  double firstStenosisLevel = ExtractStenosisLevel(data,firstDispFactor,slices,sliceAreas,data->useDiameter, data->useOldDefinition,stenosisDefs);
  double secondDispFactor = 5.0;
  double secondStenosisLevel = ExtractStenosisLevel(data,secondDispFactor,slices,sliceAreas,data->useDiameter, data->useOldDefinition,stenosisDefs);
  double currResidual = 0.0;
  double trialDispFactor = 0.0;
  double trialStenosisLevel = 0.0;
  while((!converged)&&(currIt<kMaxPicardStenosisIt)){

    // Estimate trial Diplacement Factor
    trialDispFactor = firstDispFactor + ((targetStenosisLevel-firstStenosisLevel)/(secondStenosisLevel-firstStenosisLevel))*(secondDispFactor-firstDispFactor);

    // Get trial Stenosis Level
    trialStenosisLevel = ExtractStenosisLevel(data,trialDispFactor,slices,sliceAreas,data->useDiameter, data->useOldDefinition,stenosisDefs);

    // Print slice geometry
    femUtils::PlotSlicesToVTK(std::string("slicesIt_")+boost::lexical_cast<std::string>(currIt)+std::string(".vtk"),slices);

    // Check if has converged
    currResidual = fabs(((targetStenosisLevel-trialStenosisLevel)/trialStenosisLevel));
    converged = (currResidual < kStenosisTolerance);

    // Print Iterations on screen
    printf("%5d %20e %20e %20e\n",currIt+1,trialDispFactor,trialStenosisLevel,currResidual);

    // Update Disp Factors
    firstDispFactor = secondDispFactor;
    secondDispFactor = trialDispFactor;
    // update stenosis Levels
    firstStenosisLevel = secondStenosisLevel;
    secondStenosisLevel = trialStenosisLevel;

    // Update
    currIt++;
  }
  if(!converged){
    printf("Picard Iterations not Converged!\n");
  }else{
    printf("Picard Iterations Converged!\n");
    // Print the
    if (debugMode){
      // Print slice geometry
      femUtils::PlotSlicesToVTK(std::string("modelSlices")+boost::lexical_cast<std::string>(targetStenosisLevel)+std::string(".vtk"),slices);
      // Print slice area values
      for(unsigned int loopA=0;loopA<sliceAreas.size();loopA++){
        fprintf(debugFile,"%e ",sliceAreas[loopA]);
      }
      fprintf(debugFile,"\n");
      // Print Stenosis Definitions
      fprintf(defFile,"%15.4e %15.4e %15.4e %15.4e %15.4e %15.4e\n",stenosisDefs[0],stenosisDefs[1],stenosisDefs[2],stenosisDefs[3],stenosisDefs[4],stenosisDefs[5]);

      // Close File
      fclose(debugFile);
      fclose(defFile);
    }
  }
  // Return
  return trialDispFactor;
}

// ===============================================
// EXTRACT SLICE AREA AT STENOSIS REFERENCE ORIGIN
// ===============================================
double femModel::GetReferenceArea(femInputData* data){
  // Declare
  std::vector<femModelSlice*> slices;
  // Intialize
  slices.clear();
  for(int loopA=0;loopA<3;loopA++){
    femModelSlice* slice = new femModelSlice();
    slices.push_back(slice);
  }


  // Slice Model At Current Displacement Factor
  SliceModelSkin(3,0.0,data,slices);

  // Take the area in the middle
  return slices[1]->EvalSliceArea(data->mainModelRefSystem);
}

// =====================================================
// Extract Stenosis Level for a given diplacement factor
// =====================================================
double femModel::ExtractStenosisLevel(femInputData* data, double currDispFactor,
                                      std::vector<femModelSlice*> &slices, std::vector<double> &sliceAreas,
                                      bool useDiameter, bool useOldDefinition, double* stenosisDef){

  // Declare
  std::vector<femModelSlice*> refSlices;
  double refArea = 0.0;
  double refDiam = 0.0;
  double currDiam = 0.0;
  double refNoStenosisArea = 0.0;
  double refNoStenosisDiam = 0.0;

  // Resize Vector
  // SLICES
  slices.clear();
  for(int loopA=0;loopA<kStenosisSlices;loopA++){
    femModelSlice* slice = new femModelSlice();
    slices.push_back(slice);
  }

  // Slice Model At Current Displacement Factor
  SliceModelSkin(kStenosisSlices,currDispFactor,data,slices);

  // Get Reference Area
  refArea = GetReferenceArea(data);
  refDiam = sqrt(4.0*refArea/kPI);

  // For every slice order Intersection points and Eval Area
  double maxArea = 0.0;
  double maxDiam = 0.0;
  double minArea = std::numeric_limits<double>::max();
  double minDiam = std::numeric_limits<double>::max();
  double currArea = 0.0;
  sliceAreas.clear();
  for(unsigned int loopA=0;loopA<slices.size();loopA++){
    currArea = slices[loopA]->EvalSliceArea(data->mainModelRefSystem);
    currDiam = sqrt(4.0*currArea/kPI);
    sliceAreas.push_back(currArea);
    // Store values for area
    if (currArea > maxArea){
      maxArea = currArea;
    }
    if(currArea < minArea){
      minArea = currArea;
    }
    // Store Values for Diam
    if (currDiam > maxDiam){
      maxDiam = currDiam;
    }
    if(currDiam < minDiam){
      minDiam = currDiam;
    }
  }

  // Get The Reference Area/Diameter from the Undeformed stenosis level of the Mesh
  // CAREFUL: THIS APPLIES ONLY IF THE REFERENCE STENOSISI IS SPECIFIED IN TERMS OF AREA!!!
  refNoStenosisArea = (refArea/(1.0-(data->undeformedStenosisLevel/(double)100.0)));
  refNoStenosisDiam = sqrt(4.0*refNoStenosisArea/kPI);

  // Assign Stenosis Definitions
  stenosisDef[0] = refNoStenosisDiam;
  stenosisDef[1] = refNoStenosisArea;
  stenosisDef[2] = ((1.0-(minDiam/maxDiam))*100.0);
  stenosisDef[3] = ((1.0-(minDiam/refNoStenosisDiam))*100.0);
  stenosisDef[4] = ((1.0-(minArea/maxArea))*100.0);
  stenosisDef[5] = ((1.0-(minArea/refNoStenosisArea))*100.0);

  // Eval Stenosis Level
  if(useDiameter){
    // USE DIAMETER
    if(useOldDefinition){
      return stenosisDef[2];
    }else{
      return stenosisDef[3];
    }
  }else{
    // USE AREA
    if(useOldDefinition){
      return stenosisDef[4];
    }else{
      return stenosisDef[5];
    }
  }
}

// ===============================
// Check if point is already there
// ===============================
bool pointAlreadyThere(double* backCoords,std::vector<femPoint*> &slicedPoints){
  bool found = false;
  unsigned int count = 0;
  double distance = 0.0;
  while((!found)&&(count<slicedPoints.size())){
    distance =sqrt((backCoords[0]-slicedPoints[count]->coords[0])*(backCoords[0]-slicedPoints[count]->coords[0])+
                   (backCoords[1]-slicedPoints[count]->coords[1])*(backCoords[1]-slicedPoints[count]->coords[1])+
                   (backCoords[2]-slicedPoints[count]->coords[2])*(backCoords[2]-slicedPoints[count]->coords[2]));
    // Check if found
    found = (distance < kMathZero);
    // Update
    count++;
  }
  // Return
  return found;
}

// ================
// Slice Model Skin
// ================
void femModel::SliceModelSkin(const int kStenosisSlices, double dispFactor, femInputData* data, std::vector<femModelSlice*> &slices){
  // Var
  double newCoords1[3] = {0.0};
  double newCoords2[3] = {0.0};
  double newCoords3[3] = {0.0};
  double currSliceCoord = 0.0;
  double ratio1 = 0.0;
  double ratio2 = 0.0;
  double ratio3 = 0.0;
  double backCoords[3] = {0.0};
  std::vector<femNode*> localNodeList;
  // Reset Nodes and Faces
  for(int loopC=0;loopC<kStenosisSlices;loopC++){
    slices[loopC]->slicedFaces.clear();
    slices[loopC]->slicedPoints.clear();
  }
  // Main Loop
  for(unsigned int loopA=0;loopA<faceList.size();loopA++){
    // Check if boundary Face
    if(faceList[loopA]->faceElements.size() == 1){
      // Get Face node
      nodeList[faceList[loopA]->faceNodes[0]]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,newCoords1,kDirect,kDeformed,dispFactor);
      nodeList[faceList[loopA]->faceNodes[1]]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,newCoords2,kDirect,kDeformed,dispFactor);
      nodeList[faceList[loopA]->faceNodes[2]]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,newCoords3,kDirect,kDeformed,dispFactor);
      for(int loopC=0;loopC<kStenosisSlices;loopC++){
        currSliceCoord = -0.5*data->stenosisLength + ((data->stenosisLength)/((double)(kStenosisSlices-1)))*loopC;
        ratio1 = (newCoords1[0]-currSliceCoord)/(newCoords2[0]-currSliceCoord);
        ratio2 = (newCoords2[0]-currSliceCoord)/(newCoords3[0]-currSliceCoord);
        ratio3 = (newCoords3[0]-currSliceCoord)/(newCoords1[0]-currSliceCoord);
        if((ratio1<0.0)||(ratio2<0.0)||(ratio3<0.0)){
          // Add Face
          slices[loopC]->slicedFaces.push_back(loopA);
          // Eval Intersection
          localNodeList.clear();
          femUtils::evalCoordinatePlaneIntersections(currSliceCoord,newCoords1,newCoords2,newCoords3,localNodeList);
          // Add the mapped back nodes
          for(unsigned int loopD=0;loopD<localNodeList.size();loopD++){
            // Check if it is within the stenosys Box
            if((fabs(localNodeList[loopD]->coords[1])<kStenosisBoxFactor*data->stenosisLength)&&
               (fabs(localNodeList[loopD]->coords[2])<kStenosisBoxFactor*data->stenosisLength)){
              // Transform back
              localNodeList[loopD]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,backCoords,kReverse,kUndeformed,0.0);
              // Check if point exists
              if(!pointAlreadyThere(backCoords,slices[loopC]->slicedPoints)){
                // Add to the point list
                femPoint* point = new femPoint(backCoords);
                slices[loopC]->slicedPoints.push_back(point);
              }
            }
          }
        }
      }
    }
  }
}

// ================
// WRITE CVPRE FILE
// ================
void femModel::WriteCvPreFile(std::string fileName){
  // Test Print out
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // WRITE TOTALS
  fprintf(outFile,"\n");
  fprintf(outFile,"# MODEL TOTALS\n");
  fprintf(outFile,"number_of_variables %d\n",5);
  fprintf(outFile,"number_of_nodes %d\n",(int)nodeList.size());
  fprintf(outFile,"number_of_elements %d\n",(int)elementList.size());
  fprintf(outFile,"number_of_mesh_edges %d\n",5);// !!! BOH...
  fprintf(outFile,"number_of_mesh_faces %d\n",(int)faceList.size());
  fprintf(outFile,"\n");
  fprintf(outFile,"# NODE COORDINATES AND ELEMENT CONNECTIVITY\n");
  fprintf(outFile,"nodes model.coordinates.gz\n");
  fprintf(outFile,"elements model.connectivity.gz\n");
  fprintf(outFile,"\n");
  fprintf(outFile,"# BOUNDARY FACES AND ADJIACENCY\n");
  fprintf(outFile,"boundary_faces all_exterior_elemfaces.ebc.gz\n");
  fprintf(outFile,"adjacency model.xadj.gz\n");
  fprintf(outFile,"\n");
  fprintf(outFile,"# NO SLIP\n");
  fprintf(outFile,"noslip exterior_nodes_group_0.nbc.gz\n");
  fprintf(outFile,"\n");
  fprintf(outFile,"# ZERO PRESSURE\n");
  fprintf(outFile,"zero_pressure exterior_elemfaces_group_1.ebc.gz\n");
  fprintf(outFile,"\n");
  fprintf(outFile,"# SET SURFACE ID\n");
  fprintf(outFile,"set_surface_id exterior_elemfaces_group_0.ebc.gz 0\n");
  fprintf(outFile,"set_surface_id exterior_elemfaces_group_1.ebc.gz 1\n");
  fprintf(outFile,"set_surface_id exterior_elemfaces_group_2.ebc.gz 2\n");
  fprintf(outFile,"\n");
  fprintf(outFile,"# WRITE MODEL DATA\n");
  fprintf(outFile,"write_geombc geombc.dat.1\n");
  fprintf(outFile,"write_restart restart.0.1\n");
  // Close File
  fclose(outFile);
}

// ===========================
// EXPORT FILE TO CVPRE FOLDER
// ===========================
void femModel::ExportToCvPre(double dispFactor, std::string pathName){
  // VAR
  std::string linuxCommand = "";
  int totalFaceGroups = 0;

  // CREATE MODEL FOLDER
  mkdir(pathName.c_str(),0777);

  // EXPORT COORDINATES
  std::string coordFile = pathName + std::string("/model.coordinates");
  WriteNodeCoordsToFile(dispFactor, coordFile);
  // ZIP FILE
  linuxCommand = std::string("gzip ") + coordFile;
  system(linuxCommand.c_str());

  // EXPORT CONNECTIVITIES
  std::string connFile = pathName + std::string("/model.connectivity");
  WriteElementConnectionsToFile(connFile);
  // ZIP FILE
  linuxCommand = std::string("gzip ") + connFile;
  system(linuxCommand.c_str());

  // EXPORT ADJACENCY FILE
  std::string adjFile = pathName + std::string("/model.xadj");
  WriteAdjacenciesToFile(adjFile);
  // ZIP FILE
  linuxCommand = std::string("gzip ") + adjFile;
  system(linuxCommand.c_str());

  // FORM FACE GROUPS
  FormBoundaryFaceGroups(totalFaceGroups);

  // EXPORT BOUNDARY ELEMENTS
  ExportBoundaryElementFiles(dispFactor,totalFaceGroups,pathName);
  // EXPORT BOUNDARY NODES
  ExportBoundaryNodeFiles(totalFaceGroups,pathName);

  // Write cvPre file
  std::string cvPreFile = pathName + std::string("/model.cvpre");
  WriteCvPreFile(cvPreFile);
}

// ==========================
// EXPORT ADJ FILE FOR PHASTA
// ==========================
void femModel::WriteAdjacenciesToFile(std::string fileName){
  // Test Print out
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  int adjPointer[elementList.size()+1];
  // Print total number of elements + 1
  fprintf(outFile,"xadj: %d\n",(int)elementList.size() + 1);
  int totAdj = 0;
  adjPointer[0] = totAdj;
  int currFace = 0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      if(faceList[currFace]->faceElements.size() > 1){
        totAdj++;
      }
    }
    // Update adjPointer
    adjPointer[loopA + 1] = totAdj;
  }
  // Print total adj list
  fprintf(outFile,"adjncy: %d\n",totAdj);
  // Print the Adjacency Vector
  for(unsigned int loopA=0;loopA<elementList.size() + 1;loopA++){
    fprintf(outFile,"%d\n",adjPointer[loopA]);
  }

  // Print the element Adjacencies
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      if(faceList[currFace]->faceElements.size() > 1){
        if(faceList[currFace]->faceElements[0] == (int)loopA){
          // OCCHIO!!!!! al + 1: NON CI VA?????
          fprintf(outFile,"%d\n",faceList[currFace]->faceElements[1]);
        }else{
          fprintf(outFile,"%d\n",faceList[currFace]->faceElements[0]);
        }
      }
    }
  }

  // Close File
  fclose(outFile);
}

// ====================================
// Group boundary faces based on normal
// ====================================
femModel* femModel::FormBoundaryFaceModel(){
  // Create visited node and permutation vectors
  bool isNodeVisited[(int)nodeList.size()];
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    isNodeVisited[loopA] = false;
  }
  int nodePermutation[(int)nodeList.size()];
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    nodePermutation[loopA] = -1;
  }
  // Visit the boundary faces
  int totBoundaryFaceNodes = 0;
  int currNode = 0;
  for(unsigned int loopA=0;loopA<faceList.size();loopA++){
    if(faceList[loopA]->faceElements.size() == 1){
      // loop through face nodes
      for(unsigned int loopB=0;loopB<faceList[loopA]->faceNodes.size();loopB++){
        // Current Node
        currNode = faceList[loopA]->faceNodes[loopB];
        if(!isNodeVisited[currNode]){
          // Mark as visited
          isNodeVisited[currNode] = true;
          // Store permutation
          nodePermutation[currNode] = totBoundaryFaceNodes;
          // Increment global Counter
          totBoundaryFaceNodes++;
        }
      }
    }
  }

  // Re-scan Boundary Face elements and create new model
  femModel* faceBoundaryModel = new femModel();
  femElement* newEl;
  int localNodeList[3] = {0};
  for(unsigned int loopA=0;loopA<faceList.size();loopA++){
    if(faceList[loopA]->faceElements.size() == 1){
      // Add all nodes
      for(unsigned int loopB=0;loopB<faceList[loopA]->faceNodes.size();loopB++){
        localNodeList[loopB]  = nodePermutation[faceList[loopA]->faceNodes[loopB]];
      }
      // Create new face elements
      newEl = new femTri3(loopA,1,3,localNodeList);
      // Add element to the list
      faceBoundaryModel->elementList.push_back(newEl);
    }
  }

  // Create Node List for new model
  femNode* newNode;
  double disps[6] = {0.0};
  faceBoundaryModel->nodeList.resize(totBoundaryFaceNodes);
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    if (isNodeVisited[loopA]){
      // Create New Node
      newNode = new femNode(nodePermutation[loopA],nodeList[loopA]->coords,disps);
      // Add to the list
      faceBoundaryModel->nodeList[nodePermutation[loopA]] = newNode;
    }
  }

  // Form Edge List
  faceBoundaryModel->FormElementEdgeList();

  // Return Model with only triangular faces
  return faceBoundaryModel;
}

// ================
// FORM FACE GROUPS
// ================
void femModel::FormBoundaryFaceGroups(int &totalFaceGroups){
  // Form Skin Model
  femModel* skinModel = nullptr;
  skinModel = FormBoundaryFaceModel();
  // Number Faces using Normal continuity
  skinModel->GroupFacesByNormal(totalFaceGroups);
  // Transfer face numbering
  for(unsigned int loopA=0;loopA<skinModel->elementList.size();loopA++){
    faceList[skinModel->elementList[loopA]->elementNumber]->group = skinModel->elementList[loopA]->propertyNumber;
  }
  // Export Model for debug
  skinModel->ExportToVTKLegacy("debugFaceModel.vtk");
  // Delete Temporary skin Model
  delete skinModel;
}

// ====================================================
// Export Boundary Element File for global + all groups
// ====================================================
void femModel::ExportBoundaryElementFiles(double dispFactor, int totalFaceGroups, std::string pathName){
  std::string linuxCommand;
  // Export Global Boundary Elements
  std::string currFileName = pathName + std::string("/all_exterior_elemfaces.ebc");
  ExportElementFaceGroupToFile(currFileName,-1);
  std::string vtkFileName = pathName + std::string("/all_exterior_elemfaces.vtk");
  ExportSkinFaceGroupToVTK(vtkFileName,dispFactor,-1);
  // Zip Face File
  linuxCommand = std::string("gzip ") + currFileName;
  system(linuxCommand.c_str());

  // Loop through all groups and export
  for(int loopA=0;loopA<totalFaceGroups;loopA++){
    // Export Global Boundary Elements
    currFileName = pathName + std::string("/exterior_elemfaces_group_") + boost::lexical_cast<std::string>(loopA) + std::string(".ebc");
    ExportElementFaceGroupToFile(currFileName,loopA);
    vtkFileName = pathName + std::string("/exterior_elemfaces_group_") + boost::lexical_cast<std::string>(loopA) + std::string(".vtk");
    ExportSkinFaceGroupToVTK(vtkFileName,dispFactor,loopA);
    // Zip Face File
    linuxCommand = std::string("gzip ") + currFileName;
    system(linuxCommand.c_str());
  }
}

// =================================================
// Export Boundary node File for global + all groups
// =================================================
void femModel::ExportBoundaryNodeFiles(int totalFaceGroups, std::string pathName){
  std::string linuxCommand;
  // Export Global Boundary Elements
  std::string currFileName = pathName + std::string("/all_exterior_nodes.nbc");
  ExportNodeFaceGroupToFile(currFileName,-1);
  linuxCommand = std::string("gzip ") + currFileName;
  system(linuxCommand.c_str());

  // Loop through all groups and export
  for(int loopA=0;loopA<totalFaceGroups;loopA++){
    // Export Global Boundary Elements
    currFileName = pathName + std::string("/exterior_nodes_group_") + boost::lexical_cast<std::string>(loopA) + std::string(".nbc");
    ExportNodeFaceGroupToFile(currFileName,loopA);
    linuxCommand = std::string("gzip ") + currFileName;
    system(linuxCommand.c_str());
  }
}

// =========================================
// Form Edge list for model with 2D entities
// =========================================
void femModel::FormElementEdgeList(){
  // Write Message
  femUtils::WriteMessage(std::string("Forming Edge List..."));
  int totalEdges = 0;
  int progress = 0;
  int progressCounted = 0;
  std::vector<std::vector<femEdge*>> firstNodeEdgeList;
  firstNodeEdgeList.resize((int)nodeList.size());
  std::vector<int> nodes;
  int EdgeID = -2;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    // Progress
    progress = (int)((loopA/double(elementList.size()-1))*100.0);
    if (((progress % 10) == 0)&&((progress / 10) != progressCounted)){
      progressCounted = (progress / 10);
      femUtils::WriteMessage(std::string(boost::lexical_cast<std::string>(progress))+".");
    }
    // Loop on the edges
    for(int loopB=0;loopB<elementList[loopA]->numberOfNodes;loopB++){
      nodes.clear();
      nodes.push_back(elementList[loopA]->elementConnections[loopB % (elementList[loopA]->numberOfNodes)]);
      nodes.push_back(elementList[loopA]->elementConnections[(loopB + 1) % (elementList[loopA]->numberOfNodes)]);
      // Sort Nodes
      std::sort(std::begin(nodes), std::end(nodes));
      // Check if the face is already there
      EdgeID = FindEdge(firstNodeEdgeList[nodes[0]],nodes);
      if (EdgeID>-1){
        elementList[loopA]->elementEdges.push_back(EdgeID);
        edgeList[EdgeID]->edgeElements.push_back(loopA);
      }else{
        totalEdges++;
        femEdge* edge = new femEdge();
        edgeList.push_back(edge);
        elementList[loopA]->elementEdges.push_back(totalEdges-1);
        edgeList[totalEdges-1]->edgeElements.push_back(loopA);
        for(unsigned int loopC=0;loopC<nodes.size();loopC++){
          edgeList[totalEdges-1]->edgeNodes.push_back(nodes[loopC]);
        }
        // Add to the list associated to the first node
        firstNodeEdgeList[nodes[0]].push_back(new femEdge(totalEdges-1,nodes));
      }
    }
  }
  // Write Message
  femUtils::WriteMessage(std::string("Done.\n"));
}

// ==================================
// Pick the first not visited element
// ==================================
int pickNotVisited(int size, bool* isVisited){
  int count = 0;
  bool found = false;
  while((!found)&&(count<size)){
    // Check if found
    found = (!isVisited[count]);
    // Increment counter
    if(!found){
      count++;
    }
  }
  // Return
  if(!found){
    return -1;
  }else{
    return count;
  }
}

// ========================================
// Eval element Normal only for 2D elements
// ========================================
void femModel::evalElementNormal(int firstElement, double* normal){
  if(elementList[firstElement]->is2D()){
    int elNodes[3] = {0};
    double firstVector[3] = {0.0};
    double secondVector[3] = {0.0};
    // Gather three nodes from the connectivity
    elNodes[0] = elementList[firstElement]->elementConnections[0];
    elNodes[1] = elementList[firstElement]->elementConnections[1];
    elNodes[2] = elementList[firstElement]->elementConnections[2];
    // Get the first Vector
    firstVector[0] = nodeList[elNodes[1]]->coords[0] - nodeList[elNodes[0]]->coords[0];
    firstVector[1] = nodeList[elNodes[1]]->coords[1] - nodeList[elNodes[0]]->coords[1];
    firstVector[2] = nodeList[elNodes[1]]->coords[2] - nodeList[elNodes[0]]->coords[2];
    femUtils::Normalize3DVector(firstVector);
    // Get the second Vector
    secondVector[0] = nodeList[elNodes[2]]->coords[0] - nodeList[elNodes[0]]->coords[0];
    secondVector[1] = nodeList[elNodes[2]]->coords[1] - nodeList[elNodes[0]]->coords[1];
    secondVector[2] = nodeList[elNodes[2]]->coords[2] - nodeList[elNodes[0]]->coords[2];
    femUtils::Normalize3DVector(secondVector);
    // Make External product
    femUtils::Do3DExternalProduct(firstVector,secondVector,normal);
    femUtils::Normalize3DVector(normal);
  }else{
    throw femException("Internal: Cannot Evaluate Normal of 3D element");
  }
}

// ===============================================
// Check normal compatibility between two elements
// ===============================================
bool femModel::CheckNormalCompatibility(int firstElement, int secondElement){
  double firstNormal[3] = {0.0};
  double secondNormal[3] = {0.0};
  // eval the first normal: Normalized!
  evalElementNormal(firstElement,firstNormal);
  // eval the second normal: Normalized!
  evalElementNormal(secondElement,secondNormal);

  // Get Internal Product and Check Limits
  double currProd = femUtils::Do3DInternalProduct(firstNormal,secondNormal);
  if(currProd>1.0){
    currProd = 1.0;
  }
  if(currProd<-1.0){
    currProd = -1.0;
  }
  // Get the angle in degrees
  double alpha = acos(currProd)*(180.0/kPI);
  //if (alpha > 90.0){
  //  alpha = fabs(alpha - 180.0);
  //}
  bool isCompatible = (alpha < kNormalAngleLimit);
  return isCompatible;
}

// ========================================
// GROUP FACES BY NORMAL FOR 2D SKIN MODELS
// ========================================
void femModel::GroupFacesByNormal(int &currGroup){
  // Initialize the total Number of Groups
  currGroup = 0;
  // VAR
  bool isVisited[(int)elementList.size()];
  // Initialize
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    isVisited[loopA] = false;
  }
  std::stack<int> elementStack;
  // Start from Element 0
  int currElement = 0;
  int currEdge = 0;
  int otherElement = 0;
  bool finishedAllGroups = false;
  bool finishedCurrentGroup = false;
  bool isNormalCompatible = false;
  while(!finishedAllGroups){
    // Pick the first not visited element
    currElement = pickNotVisited((int)elementList.size(),isVisited);
    if (currElement > -1){
      // Mark as visited
      isVisited[currElement] = true;
      // Assign current group number
      elementList[currElement]->propertyNumber = currGroup;
    }
    // Start Other Group
    finishedCurrentGroup = false;
    if(currElement > -1){
      while(!finishedCurrentGroup){
        // Loop through the neighbors
        for(unsigned int loopA=0;loopA<elementList[currElement]->elementEdges.size();loopA++){
          // Store current face
          currEdge = elementList[currElement]->elementEdges[loopA];
          // Get Other Element
          if(edgeList[currEdge]->edgeElements[0] == currElement){
            otherElement = edgeList[currEdge]->edgeElements[1];
          }else if(edgeList[currEdge]->edgeElements[1] == currElement){
            otherElement = edgeList[currEdge]->edgeElements[0];
          }else{
            throw femException("Internal error: Invalid Edge definition.");
          }
          // Check normal compatibility of currElement and otherElement
          isNormalCompatible = CheckNormalCompatibility(currElement,otherElement);
          if((!isVisited[otherElement])&&(isNormalCompatible)){            
            // push in stack
            elementStack.push(otherElement);
          }
        }
        // get next currElement from stack
        if(!elementStack.empty()){
          currElement = elementStack.top();
          // Mark as visited
          isVisited[currElement] = true;
          // Assign current group number
          elementList[currElement]->propertyNumber = currGroup;
          // Elimiate element at the top of the stack
          elementStack.pop();
        }else{
          // Finished with current group
          finishedCurrentGroup = true;
          // Increment group number
          currGroup++;
        }
      }
    }else{
      finishedAllGroups = true;
    }
  }
}

// ==================================
// Export Nodes for given face groups
// ==================================
void femModel::ExportNodeFaceGroupToFile(std::string fileName,int groupNumber){
  // Form node Array
  std::vector<int> bcNodeArray;
  // Assemble
  int currFace = 0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      if(((groupNumber > -1)&&(faceList[currFace]->group == groupNumber))||((groupNumber == -1)&&(faceList[currFace]->group > -1))){
        for(unsigned int loopC=0;loopC<faceList[currFace]->faceNodes.size();loopC++){
          // Add current Node
          bcNodeArray.push_back(faceList[currFace]->faceNodes[loopC]);
        }
      }
    }
  }
  // Sort Node Array
  std::sort(bcNodeArray.begin(),bcNodeArray.end());
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Print only non-multiple nodes
  int count = 0;
  int currentNode = 0;
  int nextNode = 0;
  while(count<(int)bcNodeArray.size()-1){
    // Get Nodes
    currentNode = bcNodeArray[count];
    nextNode = bcNodeArray[count+1];
    // If different then print
    if(currentNode != nextNode){
      // Print to file
      fprintf(outFile,"%d\n",currentNode + 1);
    }
    // Increment Counter
    count++;
  }
  // Print Last
  fprintf(outFile,"%d\n",bcNodeArray[count] + 1);
  // Close File
  fclose(outFile);
}

// ==================================
// Export Faces for given Face groups
// ==================================
void femModel::ExportElementFaceGroupToFile(std::string fileName,int groupNumber){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Loop through the elements
  int currFace = 0;
  int currNode = 0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    // Export All element faces having group assigned
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      if(((groupNumber > -1)&&(faceList[currFace]->group == groupNumber))||((groupNumber == -1)&&(faceList[currFace]->group > -1))){
        fprintf(outFile,"%d %d ",loopA + 1,currFace + 1);
        for(unsigned int loopC=0;loopC<faceList[currFace]->faceNodes.size();loopC++){
          // Add current Node
          currNode = faceList[currFace]->faceNodes[loopC];
          fprintf(outFile,"%d ",currNode + 1);
        }
        fprintf(outFile,"\n");
      }
    }
  }
  // Close File
  fclose(outFile);
}

// ==================================
// Check Model Minimum element Volume
// ==================================
double femModel::CheckMinimumElementVolume(double dispFactor){
  // Initialize Volume
  double minVolume = std::numeric_limits<double>::max();
  double currVolume = 0.0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    currVolume = elementList[loopA]->EvalVolume(dispFactor,nodeList);
    if(currVolume < minVolume){
      minVolume = currVolume;
    }
  }
  return minVolume;
}

// ====================================
// Check Minimum Mixed product in model
// ====================================
double femModel::CheckMinimumElementMixProduct(double dispFactor){
  double minProd = std::numeric_limits<double>::max();
  double currProd = 0.0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    currProd = elementList[loopA]->EvalMixProduct(dispFactor,nodeList);
    if(fabs(currProd) < minProd){
      minProd = fabs(currProd);
    }
  }
  return minProd;
}

// ========================
// Flip boundary face nodes
// ========================
void femModel::OrientateBoundaryFaceNodes(){
  // Scan elements
  int n0 = 0;
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  double vecA[3] = {0.0};
  double vecB[3] = {0.0};
  double vecC[3] = {0.0};
  double vec1[3] = {0.0};
  int currNode = 0;
  int currFace = 0;
  double prod = 0.0;
  int tempNode = 0;
  // Loop on elements
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(unsigned loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      // Check of boundary face
      if(faceList[currFace]->faceElements.size() == 1){
        // Get Nodes On the Face
        n0 = faceList[currFace]->faceNodes[0];
        n1 = faceList[currFace]->faceNodes[1];
        n2 = faceList[currFace]->faceNodes[2];
        // Find the last node
        bool found = false;
        int count = 0;
        while((!found)&&(count<(int)elementList[loopA]->elementConnections.size())){
          currNode = elementList[loopA]->elementConnections[count];
          found = ((currNode != n0)&&(currNode != n1)&&(currNode != n2));
          if(!found){
            // Increment counter
            count++;
          }
        }
        if(!found){
          throw femException("Internal: Could not find Node!");
        }else{
          n3 = currNode;
        }
        // Get three vectors
        vecA[0] = nodeList[n1]->coords[0] - nodeList[n0]->coords[0];
        vecA[1] = nodeList[n1]->coords[1] - nodeList[n0]->coords[1];
        vecA[2] = nodeList[n1]->coords[2] - nodeList[n0]->coords[2];
        vecB[0] = nodeList[n2]->coords[0] - nodeList[n0]->coords[0];
        vecB[1] = nodeList[n2]->coords[1] - nodeList[n0]->coords[1];
        vecB[2] = nodeList[n2]->coords[2] - nodeList[n0]->coords[2];
        vecC[0] = nodeList[n3]->coords[0] - nodeList[n0]->coords[0];
        vecC[1] = nodeList[n3]->coords[1] - nodeList[n0]->coords[1];
        vecC[2] = nodeList[n3]->coords[2] - nodeList[n0]->coords[2];
        // Get external product
        femUtils::Do3DExternalProduct(vecA,vecB,vec1);
        // Get Internal Product
        prod = femUtils::Do3DInternalProduct(vec1,vecC);
        if(prod > 0.0){
          tempNode = faceList[currFace]->faceNodes[0];
          faceList[currFace]->faceNodes[0] = faceList[currFace]->faceNodes[1];
          faceList[currFace]->faceNodes[1] = tempNode;
        }
      }
    }
  }
}

// ========================
// Export Face Group to VTK
// ========================
void femModel::ExportSkinFaceGroupToVTK(std::string fileName, double dispFactor, int groupNumber){
  std::vector<int> currFaceNodes;
  std::vector<int> currFaceElements;
  std::vector<int> compactFaceNodes;
  std::vector<int> inverseFaceNodes;
  // Assemble Temporary Vectors
  int currFace = 0;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementFaces.size();loopB++){
      currFace = elementList[loopA]->elementFaces[loopB];
      if(faceList[currFace]->faceElements.size() == 1){
        if(((groupNumber > -1)&&(faceList[currFace]->group == groupNumber))||((groupNumber == -1)&&(faceList[currFace]->group > -1))){
          // Add the three nodes in list
          for(unsigned int loopC=0;loopC<faceList[currFace]->faceNodes.size();loopC++){
            currFaceNodes.push_back(faceList[currFace]->faceNodes[loopC]);
          }
          // Add the current Face to the list
          currFaceElements.push_back(currFace);
        }
      }
    }
  }
  // Make list compact
  femUtils::MakeCompactList(currFaceNodes,compactFaceNodes);

  // Create Inverse Relationship
  femUtils::MakeInverseList(compactFaceNodes,nodeList.size(),inverseFaceNodes);

  // Open File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");

  // Write header
  fprintf(outFile,"# vtk DataFile Version 3.0\n");
  fprintf(outFile,"vtk output\n");
  fprintf(outFile,"ASCII\n");
  fprintf(outFile,"DATASET POLYDATA\n");
  fprintf(outFile,"POINTS %d float\n",(int)compactFaceNodes.size());

  // Write all node Coordinates
  int nodeCount = 0;
  for(int loopA=0;loopA<((int)compactFaceNodes.size()/3);loopA++){
    for(int loopB=0;loopB<3;loopB++){
      // Print Values
      fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[0] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[0]);
      fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[1] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[1]);
      fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[2] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[2]);
      // Update counter
      nodeCount++;
    }
    fprintf(outFile,"\n");
  }
  // Finish off
  for(int loopA=((int)compactFaceNodes.size()/3)*3;loopA<(int)compactFaceNodes.size();loopA++){
    // Print Values
    fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[0] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[0]);
    fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[1] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[1]);
    fprintf(outFile,"%e ",nodeList[compactFaceNodes[nodeCount]]->coords[2] + dispFactor * nodeList[compactFaceNodes[nodeCount]]->displacements[2]);
    // Update counter
    nodeCount++;
  }
  fprintf(outFile,"\n");

  // Build Inverse Mapping
  // OTTIMIZZARE !!!!!!!!!
  // TROPPO LUNGO!!!

  // Write element face data as polygons
  fprintf(outFile,"POLYGONS %d %d\n",(int)currFaceElements.size(),(int)currFaceElements.size()*4);
  int currNode1 = 0;
  int currNode2 = 0;
  int currNode3 = 0;
  for(unsigned int loopA=0;loopA<currFaceElements.size();loopA++){
    // Store Current Nodes
    currNode1 = faceList[currFaceElements[loopA]]->faceNodes[0];
    currNode2 = faceList[currFaceElements[loopA]]->faceNodes[1];
    currNode3 = faceList[currFaceElements[loopA]]->faceNodes[2];
    // Check if nodes are valid
    if((currNode1>-1)&&(currNode2>-1)&&(currNode3>-1)){
      fprintf(outFile,"3 %d %d %d\n",
                      inverseFaceNodes[currNode1],
                      inverseFaceNodes[currNode2],
                      inverseFaceNodes[currNode3]);
    }else{
      throw femException("Internal: Invalid Face Number.\n");
    }

  }

  // Write Node Maps to original model as scalars
  fprintf(outFile,"POINT_DATA %d\n",(int)compactFaceNodes.size());
  fprintf(outFile,"SCALARS orig_numbers int\n");
  fprintf(outFile,"LOOKUP_TABLE default\n");
  for(unsigned int loopA=0;loopA<compactFaceNodes.size();loopA++){
    fprintf(outFile,"%d \n",compactFaceNodes[loopA]);
  }

  // Close File
  fclose(outFile);
}

// =============================
// Normalize Model Displacements
// =============================
void femModel::NormalizeDisplacements(double maxDisp){
  double currMaxDisp = -std::numeric_limits<double>::max();
  double currDispModule = 0.0;
  for(unsigned loopA=0;loopA<nodeList.size();loopA++){
    // Find Max displacement in module
    currDispModule = sqrt(nodeList[loopA]->displacements[0]*nodeList[loopA]->displacements[0]+
                          nodeList[loopA]->displacements[1]*nodeList[loopA]->displacements[1]+
                          nodeList[loopA]->displacements[2]*nodeList[loopA]->displacements[2]);
    // Store Maximum Value
    if(currDispModule > currMaxDisp){
      currMaxDisp = currDispModule;
    }
  }
  // Get Scale Factor
  double scaleFactor = (double)maxDisp/(double)currMaxDisp;
  for(unsigned loopA=0;loopA<nodeList.size();loopA++){
    for(int loopB=0;loopB<6;loopB++){
      nodeList[loopA]->displacements[loopB] *= scaleFactor;
    }
  }
}

// ====================================
// EVALUATE MODEL QUALITY DISTRIBUTIONS
// ====================================
void femModel::EvalModelQualityDistributions(std::string fileName, double* limitBox){

  // Set Quantities
  int numberOfBins = 300;
  double currIntervalVolume = 0.0;
  double currIntervalMixedProduct = 0.0;
  double currVolume = 0.0;
  double currMixedProd = 0.0;
  double centroid[3] = {0.0};

  // Allocate Arrays
  double binMinVolume[numberOfBins];
  double binMaxVolume[numberOfBins];
  double binMinMixed[numberOfBins];
  double binMaxMixed[numberOfBins];
  double binCenterVolume[numberOfBins];
  double binCenterMixed[numberOfBins];
  double binArrayVolume[numberOfBins];
  double binArrayMixed[numberOfBins];

  // FORM BIN LIMITS FOR SINGLE SCAN
  // Element Volume
  femStatistics::FormBinLimits(this,kVolume,currIntervalVolume,numberOfBins,binMinVolume,binMaxVolume,binCenterVolume,limitBox);
  // Element Mixed Product
  femStatistics::FormBinLimits(this,kMixedProduct,currIntervalMixedProduct,numberOfBins,binMinMixed,binMaxMixed,binCenterMixed,limitBox);
  // Loop through all elements in the mesh
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    // Get Element Centroid
    elementList[loopA]->evalElementCentroid(nodeList,centroid);
    // Check if inside Limits
    if(femUtils::isInsideLimits(centroid,limitBox)){
      // Determine the Quality indexes
      currVolume = elementList[loopA]->EvalVolume(0.0,nodeList);
      currMixedProd = elementList[loopA]->EvalMixProduct(0.0,nodeList);
      // Assign to BIN
      // Volume
      femStatistics::AssignToBin(currVolume,numberOfBins,binMinVolume,binMaxVolume,binArrayVolume);
      // Mixed Product
      femStatistics::AssignToBin(currMixedProd,numberOfBins,binMinMixed,binMaxMixed,binArrayMixed);
    }
  }
  // NORMALIZE VALUES
  // Volume
  femStatistics::NormalizeBinArray(numberOfBins,binArrayVolume,currIntervalVolume);
  // Mixed Product
  femStatistics::NormalizeBinArray(numberOfBins,binArrayMixed,currIntervalMixedProduct);

  // Copy to Std Vector
  std::vector<double> binCenterVolumeVec;
  std::vector<double> binArrayVolumeVec;
  std::vector<double> binCenterMixedVec;
  std::vector<double> binArrayMixedVec;
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binCenterVolumeVec.push_back(binCenterVolume[loopA]);
    binArrayVolumeVec.push_back(binArrayVolume[loopA]);
    binCenterMixedVec.push_back(binCenterMixed[loopA]);
    binArrayMixedVec.push_back(binArrayMixed[loopA]);
  }

  // Write To File
  // Volume
  femUtils::WriteGraphToFile(fileName+"_VolumeDistr.dat",numberOfBins,binCenterVolumeVec,binArrayVolumeVec);
  // Mixed Product
  femUtils::WriteGraphToFile(fileName+"_MixedProdDistr.dat",numberOfBins,binCenterMixedVec,binArrayMixedVec);
}


// Check If the Limit Boxes are compatible
bool femModel::isModelCompatible(femModel* other,double tolerance){
  // Build Limit Boxes
  // Current Model
  EvalModelBox();
  // Other Models
  other->EvalModelBox();

  // Check If the Two Boxes are compatible
  return femUtils::AreCompatibleBoxes(modelBox,other->modelBox,tolerance);

}

// ================================
// READ MODEL NODES FROM VTK LEGACY
// ================================
void femModel::ReadModelNodesFromVTKFile(std::string fileName){

  // Declare input File
  std::ifstream infile;
  infile.open(fileName);

  // Declare
  std::vector<string> tokenizedString;
  int pointCount = 0;
  int nodeCount = 0;
  femNode* newNode;
  double coordX,coordY,coordZ;

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Look for the "POINTS" string
    if(tokenizedString[0] == std::string("POINTS")){
      // Read Number of Points
      pointCount = atoi(tokenizedString[1].c_str());
      for(int loopA=0;loopA<pointCount/3;loopA++){
        // Read New Line
        std::getline(infile,buffer);
        // Trim String
        boost::trim(buffer);
        // Tokenize String
        boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
        // Add Three points
        for(int loopB=0;loopB<3;loopB++){
          // Increment Node Count
          nodeCount++;
          // Constructor for femNode
          coordX = atof(tokenizedString[loopB*3].c_str());
          coordY = atof(tokenizedString[loopB*3+1].c_str());
          coordZ = atof(tokenizedString[loopB*3+2].c_str());
          newNode = new femNode(nodeCount,coordX,coordY,coordZ);
          // Add to Node List
          nodeList.push_back(newNode);
        }
      }
      // Finish Off
      if((pointCount % 3)>0){
        // Read New Line
        std::getline(infile,buffer);
        // Trim String
        boost::trim(buffer);
        // Tokenize String
        boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
        // Loop through the other ones
        for(int loopB=0;loopB<(pointCount % 3);loopB++){
          // Increment Node Count
          nodeCount++;
          // Constructor for femNode
          coordX = atof(tokenizedString[loopB*3].c_str());
          coordY = atof(tokenizedString[loopB*3+1].c_str());
          coordZ = atof(tokenizedString[loopB*3+2].c_str());
          newNode = new femNode(nodeCount,coordX,coordY,coordZ);
          // Add to Node List
          nodeList.push_back(newNode);
        }
      }
    }
  }

  // Close File
  infile.close();
}

// ===================================
// READ MODEL ELEMENTS FROM VTK LEGACY
// ===================================
void femModel::ReadModelElementsFromVTKFile(std::string fileName){

  // Declare input File
  std::ifstream infile;
  infile.open(fileName);

  // Declare
  std::vector<string> tokenizedString;
  int polygonCount = 0;
  int nodeCount = 0;
  femElement* newElement;
  int* nodeConnections;
  int currElementNumber = -1;

  // Read Data From File
  std::string buffer;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    // Look for the "POINTS" string
    if(tokenizedString[0] == std::string("POLYGONS")){
      // Read Number of Points
      polygonCount = atoi(tokenizedString[1].c_str());
      for(int loopA=0;loopA<polygonCount;loopA++){
        // Read New Line
        std::getline(infile,buffer);
        // Trim String
        boost::trim(buffer);
        // Tokenize String
        boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

        // Get The number of nodes
        nodeCount = atoi(tokenizedString[0].c_str());

        // Read Element Connections
        nodeConnections = new int[nodeCount];
        for (int loopA=0;loopA<nodeCount;loopA++){
          // node numbering starts from 1
          // nodeConnections[loopA] = atoi(tokenizedString[loopA+1].c_str()) - 1;
          nodeConnections[loopA] = atoi(tokenizedString[loopA+1].c_str());
        }

        // Create a new Element
        if(nodeCount == kTri3Nodes){
          currElementNumber++;
          newElement = new femTri3(currElementNumber-1,1,nodeCount,nodeConnections);
          // Add to Element List
          elementList.push_back(newElement);
        }else if(nodeCount == kTetra4Nodes){
          currElementNumber++;
          newElement = new femTetra4(currElementNumber-1,1,nodeCount,nodeConnections);
          // Add to Element List
          elementList.push_back(newElement);
        }else if(nodeCount == kTetra10Nodes){
          currElementNumber++;
          newElement = new femTetra10(currElementNumber-1,1,nodeCount,nodeConnections);
          // Add to Element List
          elementList.push_back(newElement);
        }else{
          throw femException("Error: Invalid Element Type");
        }
      }
    }
  }
  // Close File
  infile.close();
}

// =========================
// WRITE POLYFILE FOR TETGEN
// =========================
void femModel::WriteSkinSMeshFile(std::string polyFileName){

  // Open Output File
  FILE* outFile;
  outFile = fopen(polyFileName.c_str(),"w");

  // Part 1
  fprintf(outFile,"# Part 1 - node list\n");
  fprintf(outFile,"# node count, 3 dim, no attribute, no boundary marker\n");
  fprintf(outFile,"%d 3 0 0\n",nodeList.size());

  // Write Node Coordinates
  for(int loopA=0;loopA<nodeList.size();loopA++){
    fprintf(outFile,"%d %e %e %e\n",loopA+1,nodeList[loopA]->coords[0],nodeList[loopA]->coords[1],nodeList[loopA]->coords[2]);
  }

  // Part 2
  fprintf(outFile,"# Part 2 - facet list\n");
  fprintf(outFile,"# facet count, no boundary marker\n");
  fprintf(outFile,"%d 0\n",elementList.size());

  // Write Facets
  for(int loopA=0;loopA<elementList.size();loopA++){
    // PRINT ONLY TRI3 ELEMENTS
    if(elementList[loopA]->elementConnections.size() == kTri3Nodes){
      //fprintf(outFile,"1\n");
      fprintf(outFile,"%d %d %d %d\n",kTri3Nodes,elementList[loopA]->elementConnections[0] + 1,
                                                 elementList[loopA]->elementConnections[1] + 1,
                                                 elementList[loopA]->elementConnections[2] + 1);
    }
  }

  // Close Output file
  fclose(outFile);
}

// ================
// MESH WITH TETGEN
// ================
void femModel::MeshWithTetGen(std::string polyFileName){
  femUtils::WriteMessage(std::string("Generating Mesh with TetGen...\n"));
  std::string command(std::string("tetigen -pqQ ") + polyFileName);
  int rpl = system (command.c_str());
  if(rpl != 0){
    throw femException("Error: Cannot find tetgen.\n");
  }
}

// ===================
// Read flow rate file
// ===================
//void ReadFlowRateData(std::string flowRateFile,femFaceFlowData* data){
  // Open File
  //FILE* flowFile;
  //flowFile = fopen(flowRateFile,"w");

// Read total Groups

// Read group number

  // Close File
  //fclose(flowFile);

//}

// ================================================================
// Read time dependent flow rate file and write boundary conditions
// ================================================================
//void femModel::CreateFaceFluxFiles(std::string flowRateFile){
//  // Check if the file exists
//  bool doProcess = false;
//  FILE* flowFile;
//  flowFile = fopen("faceflow.flow","w");
//  if(flowFile != NULL){
//    doProcess = true;
//    fclose(flowFile);
//  }

//  // If file Exists then process
//  if(doProcess){
//    // Read file with flow rates
//    femFaceFlowData* data = new femTableData();
//    ReadFlowRateData(flowRateFile,data);

//    // Process all faces in the model
//    for(int loopA=0;loopA<totalFaceGroups;loopA++){
//      ProcessFaceGroupFlowRate(loopA,data[loopA]);
//    }
//  }
//}
