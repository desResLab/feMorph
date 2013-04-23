#include "femInputData.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>

femInputData::femInputData()
{
}


// Member Function
void femInputData::ReadFromFile(std::string fileName){
  // Declare input File
  std::ifstream infile;
  infile.open(fileName);

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
          // Results Mapping Model
          mappingModelResultsFileName = buffer;
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
          mainModelRefSystem[0][1] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[0][2] = atof(tokenizedString[2].c_str());
          break;
        case 8:
          // Main Model Reference System Axis 2
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelRefSystem[1][0] = atof(tokenizedString[0].c_str());
          mainModelRefSystem[1][1] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[1][2] = atof(tokenizedString[2].c_str());
          break;
        case 9:
          // Main Model Reference System Axis 3
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mainModelRefSystem[2][0] = atof(tokenizedString[0].c_str());
          mainModelRefSystem[2][1] = atof(tokenizedString[1].c_str());
          mainModelRefSystem[2][2] = atof(tokenizedString[2].c_str());
          break;
        case 10:
          // Origin on Mapping Model
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mappingModelOrigin[0] = atof(tokenizedString[0].c_str());
          mappingModelOrigin[1] = atof(tokenizedString[1].c_str());
          mappingModelOrigin[2] = atof(tokenizedString[2].c_str());
          break;
        case 11:
          // Mapping Model Reference System Axis 1
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mappingModelRefSystem[0][0] = atof(tokenizedString[0].c_str());
          mappingModelRefSystem[0][1] = atof(tokenizedString[1].c_str());
          mappingModelRefSystem[0][2] = atof(tokenizedString[2].c_str());
          break;
        case 12:
          // Mapping Model Reference System Axis 2
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mappingModelRefSystem[1][0] = atof(tokenizedString[0].c_str());
          mappingModelRefSystem[1][1] = atof(tokenizedString[1].c_str());
          mappingModelRefSystem[1][2] = atof(tokenizedString[2].c_str());
          break;
        case 13:
          // Mapping Model Reference System Axis 3
          boost::split(tokenizedString, buffer, boost::is_any_of(" ,"), boost::token_compress_on);
          mappingModelRefSystem[2][0] = atof(tokenizedString[0].c_str());
          mappingModelRefSystem[2][1] = atof(tokenizedString[1].c_str());
          mappingModelRefSystem[2][2] = atof(tokenizedString[2].c_str());
          break;
      }
    }
  }
  // Close File
  infile.close();
}
