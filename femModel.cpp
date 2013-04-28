#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "femModel.h"
#include "femGrid.h"
#include "femPoint.h"
#include "femInputData.h"
#include "femException.h"
#include "femConstants.h"
#include "femUtils.h"

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
void femModel::ReadNodeCoordsFromFile(std::string fileName){
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
  double currCoords[3] = {0.0};
  double currDisps[6] = {0.0};
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
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
void femModel::ReadElementConnectionsFromFile(std::string fileName){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Write Message
  femUtils::WriteMessage(std::string("Reading Connectivity for file ")+fileName+std::string("..."));

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  int totNodes = 0;
  int currElementNumber = 0;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
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
      throw new femException("Error: Invalid Element Type");
    }

    // Add to the node List
    elementList.push_back(newElement);
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
    fprintf(outFile,"%d ",elementList[loopA]->elementNumber);
    for(unsigned int loopB=0;loopB<elementList[loopA]->elementConnections.size();loopB++){
      fprintf(outFile,"%d ",elementList[loopA]->elementConnections[loopB]);
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
      elementID = MappingModel->FindEnclosingElementWithGrid(nodeCoords,grid->gridData[currIdx]->gridElementList);
      //femUtils::WriteMessage("IN\n");
    }else{
      elementID = -1;
      //femUtils::WriteMessage("OUTSIDE\n");
    }
    if(elementID>-1){
      // Interpolate the displacements
      MappingModel->elementList[elementID]->InterpolateElementDisplacements(nodeCoords,MappingModel->nodeList,nodeDisps);
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
    if ((newCoords[0]<(0.5*data->stenosisLength))&&(newCoords[0]>(-0.5*data->stenosisLength))){
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
    scaleFactor[loopA] = 1.2*((stenosisBox[loopA*2]-stenosisBox[loopA*2+1])/double(limRect[loopA*2]-limRect[loopA*2+1]));
  }

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
      throw new femException("Internal: Element not supported.");
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
    }else{
      fclose(outFile);
      throw new femException("Error: Invalid element to Export.");
    }
  }
  // Point Data Header
  fprintf(outFile,"POINT_DATA %d\n",int(nodeList.size()));
  // Save Displacements DX,DY,DZ as vectors
  fprintf(outFile,"VECTORS DXYZ float\n");
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
    fprintf(outFile,"%e %e %e\n",nodeList[loopA]->displacements[0],nodeList[loopA]->displacements[1],nodeList[loopA]->displacements[2]);
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
int femModel::FindEnclosingElementWithGrid(double* nodeCoords, std::vector<int> &gridElementList){
  unsigned int count = 0;
  bool found = false;
  int currElement = 0;
  while ((!found)&&(count<gridElementList.size())){
    // Store element
    currElement = gridElementList[count];
    found = elementList[currElement]->isNodeInsideElement(nodeCoords,nodeList);
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
  if (debugMode){
    debugFile = fopen("stenosisAreas.dat","a");
  }
  std::vector<femModelSlice*> slices;
  std::vector<double> sliceAreas;
  bool converged = false;
  int currIt = 0;
  double firstDispFactor = 1.0;
  double firstStenosisLevel = ExtractStenosisLevel(data,firstDispFactor,slices,sliceAreas);
  double secondDispFactor = 2.0;
  double secondStenosisLevel = ExtractStenosisLevel(data,secondDispFactor,slices,sliceAreas);
  double currResidual = 0.0;
  double trialDispFactor = 0.0;
  double trialStenosisLevel = 0.0;
  while((!converged)&&(currIt<kMaxPicardStenosisIt)){

    // Estimate trial Diplacement Factor
    trialDispFactor = firstDispFactor + ((targetStenosisLevel-firstStenosisLevel)/(secondStenosisLevel-firstStenosisLevel))*(secondDispFactor-firstDispFactor);

    // Get trial Stenosis Level
    trialStenosisLevel = ExtractStenosisLevel(data,trialDispFactor,slices,sliceAreas);

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
      fclose(debugFile);
    }
  }
  // Return
  return trialDispFactor;
}

// =====================================================
// Extract Stenosis Level for a given diplacement factor
// =====================================================
double femModel::ExtractStenosisLevel(femInputData* data, double currDispFactor, std::vector<femModelSlice*> &slices, std::vector<double> &sliceAreas){

  // Resize Vector
  slices.clear();
  for(int loopA=0;loopA<kStenosisSlices;loopA++){
    femModelSlice* slice = new femModelSlice();
    slices.push_back(slice);
  }

  // Slice Model
  SliceModelSkin(kStenosisSlices,currDispFactor,data,slices);

  // For every slice order Intersection points and Eval Area
  double maxArea = 0.0;
  double minArea = std::numeric_limits<double>::max();
  double currArea = 0.0;
  sliceAreas.clear();
  for(unsigned int loopA=0;loopA<slices.size();loopA++){
    currArea = slices[loopA]->EvalSliceArea(data->mainModelRefSystem);
    sliceAreas.push_back(currArea);
    // Store values for area
    if (currArea > maxArea){
      maxArea = currArea;
    }
    if(currArea < minArea){
      minArea = currArea;
    }
  }

  // Eval Stenosis Level
  return ((1.0-(minArea/maxArea))*100.0);
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
