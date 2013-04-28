#include <stdio.h>

#include "femModel.h"
#include "femInputData.h"
#include "femUtils.h"

int main(int argc, char **argv)
{
  // Vars
  double stenosisBox[6];

  // Set Debug Mode
  bool debugMode = true;

  // Write Application Header
  femUtils::WriteAppHeader();

  // Save Command Line Parameters and switches
  std::string paramSwitch;
  std::string paramFileName;
  if (argc == 2){
    paramSwitch = argv[1];
  } else if(argc == 3){
    paramSwitch = argv[1];
    paramFileName = argv[2];
  }else{
    // Write Help and exit
    femUtils::WriteAppHelp();
    return 0;
  }

  // Handle Switch
  if(paramSwitch == "-n"){
    femUtils::WriteMessage(std::string("Normal Execution.\n"));
    femUtils::WriteMessage(std::string("\n"));
    if(argc<3){
      femUtils::WriteMessage(std::string("Missing Input File. Please specify one.\n"));
      femUtils::WriteMessage(std::string("\n"));
      return 1;
    }
    debugMode = false;
  }else if(paramSwitch == "-d"){
    femUtils::WriteMessage(std::string("Debug Execution.\n"));
    femUtils::WriteMessage(std::string("\n"));
    if(argc<3){
      femUtils::WriteMessage(std::string("Missing Input File. Please specify one.\n"));
      femUtils::WriteMessage(std::string("\n"));
      return 1;
    }
    debugMode = true;
  }else if((paramSwitch == "-?")||(paramSwitch == "-h")){
    femUtils::WriteAppHelp();
    return 0;
  }else{
    femUtils::WriteMessage(std::string("Invalid switch. Please use switch ""-?"" for further assistance.\n"));
    femUtils::WriteMessage(std::string("\n"));
    return 1;
  }


  // Read All input parameters
  femInputData* data = new femInputData();
  data->ReadFromFile(paramFileName);

  // Export Main File for Debug
  if (debugMode){
    femUtils::ExportReferenceToVTKLegacy(data,std::string("mainreference.vtk"));
  }

  // Read Main Model
  femModel* mainModel = new femModel();
  mainModel->ReadNodeCoordsFromFile(data->mainModelCoordsFileName);
  mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName);
  // Form Face List
  mainModel->FormElementFaceList();

  // Export Main File for Debug
  if (debugMode){
    mainModel->ExportToVTKLegacy(std::string("mainmodel.vtk"));
  }

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName);
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName);
  // Form Face List
  mappingModel->FormElementFaceList();

  // Read Displacement Without Rotations
  mappingModel->ReadNodeDisplacementsFromFile(data->mappingModelResultsFileName, false);

  // Export Mapping Model for Debug
  if (debugMode){
    mappingModel->ExportToVTKLegacy(std::string("mappingmodel.vtk"));
  }

  // Get Stenosis Box from main model
  mainModel->GetStenosisBox(data,stenosisBox);
  std::vector<femNode*> steNodeList;
  mainModel->createStenosisNodeList(stenosisBox,data,steNodeList);
  // Get The Centre of rotate Stenosis Box
  double stenosisBoxCenter[3] = {0.0};
  femUtils::GetAverageNodeFromList(steNodeList,stenosisBoxCenter);

  // Export Mapping Model for Debug
  if (debugMode){
    femUtils::ExportStenosisBoxToVTK(std::string("stenosisBox.vtk"),steNodeList);
  }

  // Transform Mapping Model and Displacements
  femModel* newModel;
  newModel = mappingModel->TransformModel(data,stenosisBoxCenter,stenosisBox);

  // Export Mapping Model for Debug
  if (debugMode){
    newModel->ExportToVTKLegacy(std::string("transformedmappingmodel.vtk"));
  }

  // Map Models
  double dispScaleFactor = 1.0;
  mainModel->MapDisplacements(newModel, data, dispScaleFactor);

  // Export Mapping Model for Debug
  if (debugMode){
    mainModel->ExportToVTKLegacy(std::string("finalmodel.vtk"));
  }

  // Create models with given stenotic target
  // and write them to separate files
  // Get Corresponding Displacement Scaling Factor
  FILE* areaFile;
  if(debugMode){
    areaFile = fopen("stenosisAreas.dat","w");
    fclose(areaFile);
  }
  double currDispFactor = 1.0;
  double currStenosisLevel = 0.0;
  for(unsigned int loopA=0;loopA<data->stenosisLevels.size();loopA++){
    // Get current Stenotic level
    currStenosisLevel = data->stenosisLevels[loopA];

    // Eval the stenotic displacement factor
    currDispFactor = mainModel->seekStenoticDisplacementFactor(data,currStenosisLevel,debugMode);

    // Write Output Model
    // Change Name
    mainModel->WriteNodeCoordsToFile(currDispFactor,data->mainModelCoordsOutputName);
    mainModel->WriteElementConnectionsToFile(data->mainModelConnectionsOutputName);
  }

  // Write Application Header
  femUtils::WriteMessage(std::string("\n"));
  femUtils::WriteMessage(std::string("Completed!\n"));

  // Deallocate Pointers
  delete data;
  delete mainModel;
  delete mappingModel;
  delete newModel;

  // Done
  return 0;
}
