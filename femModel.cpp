#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "femModel.h"
#include "femInputData.h"
#include "femException.h"
#include "femConstants.h"
#include "femUtils.h"

using namespace std;

// Constructor
femModel::femModel()
{
  // Initialize Model Box
  double* modelBox = new double[6];
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
}

// Distructor
femModel::~femModel()
{
  delete [] modelBox;
}

// ====================================================
// Eval Distance between the centroid and a given point
// ====================================================
double femModel::evalPointToElementDistance(int elementID, double* pointCoords){
  double centroid[3] = {0.0};
  evalElementCentroid(elementID,centroid);
  double dist = 0.0;
  for(int loopA=0;loopA<3;loopA++){
    dist += (centroid[loopA]-pointCoords[loopA])*(centroid[loopA]-pointCoords[loopA]);
  }
  // Return
  return sqrt(dist);
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
    femUtils::Rotate3DVectorAroundAxis(nodeList[loopA]->coords,angle,axis);
    // Make local Copy
    for(int loopA=0;loopA<3;loopA++){
      dispVec[loopA] = nodeList[loopA]->displacements[loopA];
      rotVec[loopA] = nodeList[loopA]->displacements[loopA+3];
    }
    // Rotate Displacements and Rotations
    femUtils::Rotate3DVectorAroundAxis(dispVec,angle,axis);
    femUtils::Rotate3DVectorAroundAxis(rotVec,angle,axis);
    // Copy back
    for(int loopA=0;loopA<3;loopA++){
      nodeList[loopA]->displacements[loopA] = dispVec[loopA];
      nodeList[loopA]->displacements[loopA+3] = rotVec[loopA];
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

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  // Initialize
  int currNodeNumber = 0;
  double currNodeX = 0.0;
  double currNodeY = 0.0;
  double currNodeZ = 0.0;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

    // Read node Number
    currNodeNumber = atoi(tokenizedString[1].c_str());

    // Read Coordinates
    currNodeX = atof(tokenizedString[2].c_str());
    currNodeY = atof(tokenizedString[3].c_str());
    currNodeZ = atof(tokenizedString[4].c_str());

    // Create a new Node
    femNode* newNode = new femNode(currNodeNumber,currNodeX,currNodeY,currNodeZ);

    // Add to the node List
    nodeList.push_back(newNode);
  }

  // Close File
  infile.close();
}

// =====================================
// Read Element Connectivities From File
// =====================================
void femModel::ReadElementConnectionsFromFile(std::string fileName){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

  // Read Data From File
  std::string buffer;
  std::vector<string> tokenizedString;
  int currElementNumber = 0;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Tokenize String
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);

    // Read node Number
    currElementNumber = atoi(tokenizedString[1].c_str());

    // Read Element Connections
    int* nodeConnections = new int(tokenizedString.size()-1);
    for (unsigned int loopA=0;loopA<tokenizedString.size()-1;loopA++){
      nodeConnections[loopA] = atoi(tokenizedString[loopA+1].c_str());
    }

    // Create a new Element: Fixed Property Number for the time being!!!
    femElement* newElement = new femElement(currElementNumber,1,tokenizedString.size()-1,nodeConnections);

    // Add to the node List
    elementList.push_back(newElement);
  }

  // Close File
  infile.close();
}

// =================================
// Read Node Displacements From File
// =================================
void femModel::ReadNodeDisplacementsFromFile(std::string fileName, bool readRotations){
  // Declare input File
  ifstream infile;
  infile.open(fileName);

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
    currNodeNumber = atoi(tokenizedString[1].c_str());

    // Read Coordinates
    currNodeDX = atof(tokenizedString[2].c_str());
    currNodeDY = atof(tokenizedString[3].c_str());
    currNodeDZ = atof(tokenizedString[4].c_str());

    // Read rotations if required
    if (readRotations) {
      currNodeRX = atof(tokenizedString[5].c_str());
      currNodeRY = atof(tokenizedString[6].c_str());
      currNodeRZ = atof(tokenizedString[7].c_str());
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
}

// ====================
// Write Coords To File
// ====================
void femModel::WriteNodeCoordsToFile(std::string fileName){
  // Open Output File
  FILE* outFile;
  outFile = fopen(fileName.c_str(),"w");
  // Write Header
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){
      fprintf(outFile,"%d %16.8e %16.8e %16.8e\n",nodeList[loopA]->nodeNumber,
              nodeList[loopA]->coords[0],nodeList[loopA]->coords[1],nodeList[loopA]->coords[2]);
  }
  // Close Output file
  fclose(outFile);
}

// =================================
// Write Element Connections To File
// =================================
void femModel::WriteElementConnectionsToFile(std::string fileName){
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
}

// ===================================================
// Map Displacements between a main and mapping models
// ===================================================
void femModel::MapDisplacements(femModel* MappingModel,
                                femInputData* data,
                                double dispScaleFactor){
  // Initialize Mapping Coords
  double NodeMappingDisplacements[3] = {0.0};
  double nodeDisps[3] = {0.0};
  // Loop through the nodes in the main model
  int elementID = 0;
  for(unsigned int loopA=0;loopA<MappingModel->nodeList.size();loopA++){
    // Find the enclosing element
    elementID = MappingModel->FindEnclosingElement(nodeList[loopA]->coords);
    // Interpolate the displacements
    elementList[elementID]->InterpolateElementDisplacements(nodeList[loopA]->coords,nodeList,NodeMappingDisplacements);
    // Save Displacements
    nodeList[loopA]->setDisplacements(nodeDisps[0],nodeDisps[1],nodeDisps[2],0.0,0.0,0.0);
  }
}

// ==================================
// Find Face Given the Attached Nodes
// ==================================
int FindFace(std::vector<femFace*> &faceList,std::vector<int> nodes){
  // Sort Nodes
  std::sort(std::begin(nodes), std::end(nodes));
  unsigned int count = 0;
  bool hasSameNodes = false;
  bool found = false;
  while ((!found)&&(count<faceList.size())){
    hasSameNodes = false;
    for(unsigned int loopA=0;loopA<faceList[count]->faceNodes.size();loopA++){
      hasSameNodes = (hasSameNodes)&&(faceList[count]->faceNodes[loopA] == nodes[loopA]);
    }
    // Assign Found
    found = hasSameNodes;
    // Update
    count++;
  }
  if(found){
    return (count-1);
  }else{
    return -1;
  }
}

// ======================
// Form Element-Face List
// ======================
// WARNING: WORKS ONLY FOR TETRA4 AND TETRA10!!!
void femModel::FormElementFaceList(){
  int totalFaces = 0;
  std::vector<int> nodes;
  int faceID = -2;
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    for(int loopB=0;loopB<kTetraFaces;loopB++){
      nodes.push_back(elementList[loopA]->elementConnections[loopB % (kTetraFaces-1)]);
      nodes.push_back(elementList[loopA]->elementConnections[loopB % (kTetraFaces-1)]);
      nodes.push_back(elementList[loopA]->elementConnections[loopB % (kTetraFaces-1)]);
      // Check if the face is already there
      faceID = FindFace(faceList,nodes);
      if (faceID>-1){
        elementList[loopA]->elementFaces.push_back(faceID);
        faceList[faceID]->faceElements.push_back(loopA);
        for(unsigned int loopC=0;loopC<nodes.size();loopC++){
          faceList[faceID]->faceNodes.push_back(nodes[loopC]);
        }
      }else{
        totalFaces++;
        femFace* face = new femFace(nodes);
        faceList.push_back(face);
        elementList[loopA]->elementFaces.push_back(totalFaces);
        faceList[totalFaces]->faceElements.push_back(loopA);
        for(unsigned int loopC=0;loopC<nodes.size();loopC++){
          faceList[totalFaces]->faceNodes.push_back(nodes[loopC]);
        }
      }
    }
  }
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
    for(int loopB=0;loopB<2;loopB++){
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
    currDist = evalPointToElementDistance(searchedElements[loopA],nodeCoords);
    if(currDist<nextDistance){
      nextDistance = currDist;
      nextElement = searchedElements[loopA];
    }
  }
}

// ==========================
// Find the Enclosing element
// ==========================
int femModel::FindEnclosingElement(double* nodeCoords){
  // Find the Enclosing element for the node
  bool isOutside = IsOutsideLimits(nodeCoords);
  if(isOutside){
    return -1;
  }else{
    // Initialize
    bool found = false;
    int countMC = 0;
    int currElement = 0;
    double currDistance = 0.0;
    int nextElement = 0;
    double nextDistance = 0.0;
    bool endSearch = false;
    // Loop on a Monte Carlo possible choice of nodes
    while ((!found)&&(countMC<kMaxEnclosingMC)){
      // Start from a random Element
      currElement = femUtils::GenerateUniformIntegers(0,elementList.size());
      currDistance = evalPointToElementDistance(currElement,nodeCoords);
      // Search for a local miminum in the distance function
      endSearch = false;
      while(!endSearch){
        // Get Next Element
        getNextElement(currElement,nodeCoords,nextElement,nextDistance);
        endSearch = (nextDistance>currDistance);
        // Update
        currDistance = nextDistance;
        currElement = nextElement;
      }

      // Check if the Element has been found
      found = elementList[currElement]->isNodeInsideElement(nodeCoords);

      // Update
      countMC++;
    }
    // Return from function
    return currElement;
  }
}

// =====================
// Eval Element Centroid
// =====================
void femModel::evalElementCentroid(int elementID, double* centroid){
  // Initialize Centroid
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  int currNode = 0;
  int totNodes = elementList[elementID]->elementConnections.size();
  for(int loopA=0;loopA<totNodes;loopA++){
    currNode = elementList[elementID]->elementConnections[loopA];
    centroid[0] = centroid[0] + nodeList[currNode]->coords[0];
    centroid[1] = centroid[1] + nodeList[currNode]->coords[1];
    centroid[2] = centroid[2] + nodeList[currNode]->coords[2];
  }
  centroid[0] = centroid[0]/double(totNodes);
  centroid[1] = centroid[1]/double(totNodes);
  centroid[2] = centroid[2]/double(totNodes);
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
bool femModel::IsOutsideLimits(double* nodeCoords){
  bool isInside = false;
  isInside = ((nodeCoords[0]>=modelBox[0])&&(nodeCoords[0]<=modelBox[0])&&
             (nodeCoords[1]>=modelBox[1])&&(nodeCoords[1]<=modelBox[1])&&
             (nodeCoords[2]>=modelBox[2])&&(nodeCoords[2]<=modelBox[2]));
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
  limRect[0] = std::numeric_limits<double>::max();
  limRect[1] = -std::numeric_limits<double>::max();
  limRect[2] = std::numeric_limits<double>::max();
  limRect[3] = -std::numeric_limits<double>::max();
  limRect[4] = std::numeric_limits<double>::max();
  limRect[5] = -std::numeric_limits<double>::max();

  // Loop through all the nodes
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){

    // Transform Node Coords
    nodeList[loopA]->TransformNodeCoords(data->mainModelOrigin,data->mainModelRefSystem,newCoords);

    // Set limits in the X direction
    limRect[0] = -0.5 * data->stenosisLength;
    limRect[1] =  0.5 * data->stenosisLength;

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

// Store Limits
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
  angle2 = femUtils::Do3DInternalProduct(mapVec2,secondVec);
}

// =========================================
// Transform node position and displacements
// =========================================
femModel* femModel::TransformModel(femInputData* data, double* stenosisBox){

  // Create the new Model
  femModel* target = new femModel();

  // Copy everything but the nodes from the Current Model
  CopyElementsTo(target);
  CopyFacesTo(target);
  CopyPropertyTo(target);

  // Intialize Limits
  double limRect[6];
  limRect[0] = std::numeric_limits<double>::max();
  limRect[1] = -std::numeric_limits<double>::max();
  limRect[2] = std::numeric_limits<double>::max();
  limRect[3] = -std::numeric_limits<double>::max();
  limRect[4] = std::numeric_limits<double>::max();
  limRect[5] = -std::numeric_limits<double>::max();

  // Transform in its own system
  double newCoords[3];
  for(unsigned int loopA=0;loopA<nodeList.size();loopA++){

    // Transform Node Coords
    nodeList[loopA]->TransformNodeCoords(data->mappingModelOrigin,data->mappingModelRefSystem,newCoords);

    // Transform Displacements
    double newDisps[3];
    nodeList[loopA]->TransformDisplacements(data->mainModelRefSystem,newDisps);

    // Add to the Target Node List
    femNode* node = new femNode(loopA,newCoords,newDisps);
    target->nodeList.push_back(node);

    // Store Limits
    updateLimits(newCoords,limRect);
  }

  // Find Scale Factors
  double scaleFactor[3];
  for(int loopA=0;loopA<3;loopA++){
    scaleFactor[loopA] = (stenosisBox[loopA*2]-stenosisBox[loopA*2+1])/double(limRect[loopA*2]-limRect[loopA*2+1]);
  }

  // Stretch Model: (0,0,0) is the new Origin!!!!
  // Also move it to the new location
  for(unsigned int loopA=0;loopA<target->nodeList.size();loopA++){
    for(int loopB=0;loopB<3;loopB++){
      // Scale and translate it to the new location
      target->nodeList[loopA]->coords[loopB] = data->mainModelOrigin[loopB] + target->nodeList[loopA]->coords[loopB]*scaleFactor[loopB];
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

  // Return Transformed Model
  return target;
}

// Copy elements to other model
void femModel::CopyElementsTo(femModel* otherModel){
  for(unsigned int loopA=0;loopA<elementList.size();loopA++){
    femElement* newElement = new femElement(elementList[loopA]);
    otherModel->elementList.push_back(newElement);
  }
}

// Copy faces to other model
void femModel::CopyFacesTo(femModel* otherModel){
  for(unsigned int loopA=0;loopA<faceList.size();loopA++){
    femFace* newFace = new femFace(faceList[loopA]);
    otherModel->faceList.push_back(newFace);
  }
}

// Copy Property to other model
void femModel::CopyPropertyTo(femModel* otherModel){
  for(unsigned int loopA=0;loopA<propList.size();loopA++){
    femProperty* newProp = new femProperty(propList[loopA]);
    otherModel->propList.push_back(newProp);
  }
}
