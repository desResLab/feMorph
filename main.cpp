#include <stdio.h>

#include "femModel.h"
#include "femInputData.h"

int main(int argc, char **argv)
{
  // Vars
  double stenosisBox[6];

  // Save Command Line Parameters as String
  std::string paramFileName(argv[1]);

  // Read All input parameters
  femInputData* data = new femInputData();
  data->ReadFromFile(paramFileName);

  // Read Main Model
  femModel* mainModel = new femModel();
  mainModel->ReadNodeCoordsFromFile(data->mainModelCoordsFileName);
  mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName);

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName);
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName);
  // Read Displacement Without Rotations
  mappingModel->ReadNodeDisplacementsFromFile(data->mappingModelResultsFileName, false);

  // Get Stenosis Box from main model
  mainModel->GetStenosisBox(data,stenosisBox);

  // Transform Mapping Model and Displacements
  femModel* newModel;
  newModel = mappingModel->TransformModel(data,stenosisBox);

  // Map Models
  double dispScaleFactor = 1.0;
  mainModel->MapDisplacements(newModel, data, dispScaleFactor);

  // Write Output Model
  mainModel->WriteNodeCoordsToFile(data->mainModelCoordsOutputName);
  mainModel->WriteElementConnectionsToFile(data->mainModelConnectionsOutputName);

  // Done
  return 0;
}
