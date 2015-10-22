#ifndef FEMINPUTDATA_H
#define FEMINPUTDATA_H

#include <stdio.h>
#include <string>
#include <vector>

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
  // =====
  // OTHER
  // =====
  // Reference Length
  double stenosisLength;
  double undeformedStenosisLevel;
  std::vector<double> stenosisLevels;
  // Type of Stenosis
  bool useDiameter;
  bool useOldDefinition;
  // Radius of circular stenosis
  double mappingDisplacementRadius;
  bool accountForAngle;
  double maxAngle;

  // Constructor and Destructor
  femInputData();
  ~femInputData();
  // Member Function
  void ReadFromFile(std::string fileName);
  // Check Normalization
  void checkAxisNormalization();
};

#endif // FEMINPUTDATA_H
