#ifndef FEMINPUTDATA_H
#define FEMINPUTDATA_H

#include <stdio.h>
#include <string>

class femInputData
{
public:
  // Data Members
  // =====
  // FILES
  // =====
  // File Names: Main Model
  std::string mainModelCoordsFileName;
  std::string mainModelConnectionsFileName;
  // File Names: Mapping Model
  std::string mappingModelCoordsFileName;
  std::string mappingModeConnectionslFileName;
  std::string mappingModelResultsFileName;
  // Output File Names
  std::string mainModelCoordsOutputName;
  std::string mainModelConnectionsOutputName;
  // ================
  // REFERENCE SYSTEM
  // ================
  // Orientation on mail model
  double* mainModelOrigin;
  double** mainModelRefSystem;
  // Orientation on mapping model
  double* mappingModelOrigin;
  double** mappingModelRefSystem;
  // =====
  // OTHER
  // =====
  // Reference Length
  double stenosisLength;

  // Constructor and Destructor
  femInputData();
  ~femInputData();
  // Member Function
  void ReadFromFile(std::string fileName);
};

#endif // FEMINPUTDATA_H
