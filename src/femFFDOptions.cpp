#include "femFFDOptions.h"

femFFDOptions::femFFDOptions(){
}

femFFDOptions::~femFFDOptions(){
}

// Read from file
void femFFDOptions::readFromFile(string fileName){
  // Write Message
  printf("Reading FFD Input file: %s\n",fileName.c_str());

  // Declare input File
  std::ifstream infile;
  infile.open(fileName.c_str());

  // Declare
  vector<string> tokenizedString;
  femDoubleVec temp;

  // Read Data From File
  std::string buffer;

  // Read the VTK input file
  std::getline(infile,buffer);
  boost::trim(buffer);
  inputFileName = buffer;


  // Read The Number of PointGrids
  std::getline(infile,buffer);
  boost::trim(buffer);
  int totGrids = atoi(buffer.c_str());

  // Loop on the total Grids
  ffdGrid locGrid;
  int totDisps = 0;
  for(int loopA=0;loopA<totGrids;loopA++){
    // Read Origin Point
    std::getline(infile,buffer);
    boost::trim(buffer);
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    locGrid.minPoint[0] = atof(tokenizedString[0].c_str());
    locGrid.minPoint[1] = atof(tokenizedString[1].c_str());
    locGrid.minPoint[2] = atof(tokenizedString[2].c_str());
    // Read Total Number of Points
    std::getline(infile,buffer);
    boost::trim(buffer);
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    locGrid.totPoints[0] = atoi(tokenizedString[0].c_str());
    locGrid.totPoints[1] = atoi(tokenizedString[1].c_str());
    locGrid.totPoints[2] = atoi(tokenizedString[2].c_str());
    // Read the FFD axis system
    std::getline(infile,buffer);
    boost::trim(buffer);
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    locGrid.gridAxis_S[0] = atof(tokenizedString[0].c_str());
    locGrid.gridAxis_S[1] = atof(tokenizedString[1].c_str());
    locGrid.gridAxis_S[2] = atof(tokenizedString[2].c_str());
    std::getline(infile,buffer);
    boost::trim(buffer);
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    locGrid.gridAxis_T[0] = atof(tokenizedString[0].c_str());
    locGrid.gridAxis_T[1] = atof(tokenizedString[1].c_str());
    locGrid.gridAxis_T[2] = atof(tokenizedString[2].c_str());
    std::getline(infile,buffer);
    boost::trim(buffer);
    boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
    locGrid.gridAxis_U[0] = atof(tokenizedString[0].c_str());
    locGrid.gridAxis_U[1] = atof(tokenizedString[1].c_str());
    locGrid.gridAxis_U[2] = atof(tokenizedString[2].c_str());
    // Read the total number of displacements
    std::getline(infile,buffer);
    boost::trim(buffer);
    totDisps = atoi(buffer.c_str());
    locGrid.dispNodes.clear();
    locGrid.dispVals.clear();
    // Read the displacements
    for(int loopB=0;loopB<totDisps;loopB++){
      std::getline(infile,buffer);
      boost::trim(buffer);
      boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
      locGrid.dispNodes.push_back(atoi(tokenizedString[0].c_str()));
      temp.clear();
      temp.push_back(atof(tokenizedString[1].c_str()));
      temp.push_back(atof(tokenizedString[2].c_str()));
      temp.push_back(atof(tokenizedString[3].c_str()));
      locGrid.dispVals.push_back(temp);
    }
    // Add to the vector
    ffdData.push_back(locGrid);
  }
  // Close File
  infile.close();
}

