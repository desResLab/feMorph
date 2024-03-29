#include "femInputData.h"
#include "femUtils.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>

// Contructor
femInputData::femInputData()
{
  // Allocate
  mainModelOrigin = new double[3];
  mainModelRefSystem = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    mainModelRefSystem[loopA] = new double[3];
  }
  stenosisLength = 0.0;
}

// Destructor
femInputData::~femInputData(){
  // DeAllocate
  delete [] mainModelOrigin;
  for(int loopA=0;loopA<3;loopA++){
    delete [] mainModelRefSystem[loopA];
  }
  delete [] mainModelRefSystem;
}


// Member Function
void femInputData::ReadFromFile(std::string fileName){
  // Declare input File
  std::ifstream infile;
  infile.open(fileName.c_str());

  // Write Message
  femUtils::WriteMessage("Reading Input Parameter...\n");

  // Read Data From File
  int lineCount = 0;
  std::string buffer;
  std::vector<std::string> tokenizedString;
  while (std::getline(infile,buffer)){
    // Trim String
    boost::trim(buffer);
    // Avoid Comments On File
    if (buffer[0] != '#'){
      // Update Line Count
      lineCount++;
      // Choose What to Read
      switch (lineCount){
        case 1:
          // Coordinates Main Model
          mainModelCoordsFileName = buffer;
          break;
        case 2:
          // Connections Main Model
          mainModelConnectionsFileName = buffer;
          break;
        case 3:
          // Coordinates Mapping Model
          mappingModelCoordsFileName = buffer;
          break;
        case 4:
          // Connections Mapping Model
          mappingModeConnectionslFileName = buffer;
          break;
        case 5:
          // Displacements of the mapping model
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          if(boost::to_upper_copy(tokenizedString[0]) == "USEFILE"){
            mappingDisplacementType = ipUseFile;
            mappingModelResultsFileName = tokenizedString[1];
          }else if(boost::to_upper_copy(tokenizedString[0]) == "USEPARAMS"){
            mappingDisplacementType = ipUseParams;
            mappingDisplacementRadius = atof(tokenizedString[1].c_str());
            if(tokenizedString[2] == "1"){
              accountForAngle = true;
              maxAngle = atof(tokenizedString[3].c_str());
            }else{
              accountForAngle = false;
            }
          }
          break;
        case 6:
          // Origin On Main Model
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelOrigin[0] = atof(tokenizedString[0].c_str());
          mainModelOrigin[1] = atof(tokenizedString[1].c_str());
          mainModelOrigin[2] = atof(tokenizedString[2].c_str());
          break;
        case 7:
          // Main Model Reference System Axis 1
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelRefSystem[0][0] = atof(tokenizedString[0].c_str());
          mainModelRefSystem[1][0] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[2][0] = atof(tokenizedString[2].c_str());
          break;
        case 8:
          // Main Model Reference System Axis 2
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelRefSystem[0][1] = atof(tokenizedString[0].c_str());
          mainModelRefSystem[1][1] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[2][1] = atof(tokenizedString[2].c_str());
          break;
        case 9:
          // Main Model Reference System Axis 3
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelRefSystem[0][2] = atof(tokenizedString[0].c_str());
          mainModelRefSystem[1][2] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[2][2] = atof(tokenizedString[2].c_str());
          break;
        case 10:
          // Length of the stenosis
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          stenosisLength = atof(tokenizedString[0].c_str());
          break;
        case 11:
          // Levels of the stenosis
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          for(unsigned int loopA=0;loopA<tokenizedString.size();loopA++){
            stenosisLevels.push_back(atof(tokenizedString[loopA].c_str()));
          }
          break;
        case 12:
          // Undeformed Stenosis level, use diameter and use Old definition
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          undeformedStenosisLevel = atof(tokenizedString[0].c_str());
          // Type of Stenosis
          useDiameter= (atoi(tokenizedString[1].c_str()) == 1);
          useOldDefinition = (atoi(tokenizedString[2].c_str()) == 1);
          break;
      }
    }
  }
  // Set the name of the output file
  mainModelCoordsOutputName = std::string("result_model.coordinates");
  mainModelConnectionsOutputName = std::string("result_model.connections");

  // Check Axis Normalization
  checkAxisNormalization();

  // Close File
  infile.close();

  femUtils::WriteMessage("Done.\n");
}

// =========================================
// RENORMALIZE AXIS SYSTEM AND CHECK RH RULE
// =========================================
void femInputData::checkAxisNormalization(){
  double currMod = 0.0;
  for(int loopA=0;loopA<kDims;loopA++){
    currMod = 0.0;
    for(int loopB=0;loopB<kDims;loopB++){
      currMod += mainModelRefSystem[loopB][loopA]*mainModelRefSystem[loopB][loopA];
    }
    currMod = sqrt(currMod);
    if(abs(currMod-1.0)>kMathZero){
      femUtils::WriteMessage(std::string("Warning: Axis ") + femUtils::intToStr(loopA+1) + std::string(" was normalized.\n"));
      for(int loopB=0;loopB<kDims;loopB++){
        mainModelRefSystem[loopB][loopA] = mainModelRefSystem[loopB][loopA]/currMod;
      }
    }
  }

  // Check RH Rule
  double axis1[3] = {0.0};
  double axis2[3] = {0.0};
  double axis3[3] = {0.0};
  for(int loopA=0;loopA<kDims;loopA++){
    axis1[loopA] = mainModelRefSystem[loopA][0];
    axis2[loopA] = mainModelRefSystem[loopA][1];
  }
  femUtils::Do3DExternalProduct(axis1,axis2,axis3);
  double intProd = axis3[0]*mainModelRefSystem[0][2] +
                   axis3[1]*mainModelRefSystem[1][2] +
                   axis3[2]*mainModelRefSystem[2][2];
  if(intProd<0.0){
    throw femException("The axis system is not consistent with the right-hand rule.\n");
  }

}
