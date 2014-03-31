#include <stdio.h>
#include <boost/lexical_cast.hpp>

#include "femModel.h"
#include "femInputData.h"
#include "femUtils.h"
#include "femWeightedFileName.h"
#include "femModelSequence.h"
#include "femProgramOptions.h"

// ====================
// Normal Model Running
// ====================
int runNormalMode(femProgramOptions* options){
  // Vars
  double stenosisBox[6];

  // Set Debug Mode
  bool debugMode = true;
  bool reducedOutput = true;

  // Read All input parameters
  femInputData* data = new femInputData();
  data->ReadFromFile(options->inputFileName);

  // Export Main File for Debug
  if (debugMode){
    femUtils::ExportReferenceToVTKLegacy(data,std::string("mainreference.vtk"));
  }

  // Read Main Model
  femModel* mainModel = new femModel();
  mainModel->ReadNodeCoordsFromFile(data->mainModelCoordsFileName,false);
  mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName,false);
  // Form Face List
  mainModel->FormElementFaceList();

  // Export Main File for Debug
  if (debugMode){
    mainModel->ExportToVTKLegacy(std::string("mainmodel.vtk"));
  }

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName,false);
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName,false);
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

  // Normalize Model Displacements
  mainModel->NormalizeDisplacements(0.01*data->stenosisLength);

  // Export Mapping Model for Debug
  if (debugMode){
    mainModel->ExportToVTKLegacy(std::string("finalmodel.vtk"));
  }

  // Orientate face nodes before exporting
  mainModel->OrientateBoundaryFaceNodes();

  // Export cvPre Model ready to presolve
  if(!reducedOutput){
    mainModel->ExportToCvPre(0.0,std::string("reference"));
  }

  // Check Volume of Undeformed Mesh
  double minVol = mainModel->CheckMinimumElementVolume(0.0);
  femUtils::WriteMessage(std::string("Undeformed - Minimum element Volume: ") + boost::lexical_cast<std::string>(minVol) + std::string("\n"));
  double minMixedProduct = mainModel->CheckMinimumElementMixProduct(0.0);
  femUtils::WriteMessage(std::string("Undeformed - Minimum element Mixed Product: ") + boost::lexical_cast<std::string>(minMixedProduct) + std::string("\n"));

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

    // Write Message
    femUtils::WriteMessage(std::string("Computing Stenosis Level ") + boost::lexical_cast<std::string>(currStenosisLevel) + std::string("...\n"));

    // Eval the stenotic displacement factor
    currDispFactor = mainModel->seekStenoticDisplacementFactor(data,currStenosisLevel,debugMode);

    // Check Volume of Undeformed Mesh
    double minVol = mainModel->CheckMinimumElementVolume(currDispFactor);
    femUtils::WriteMessage(std::string("Minimum element Volume: ") + boost::lexical_cast<std::string>(minVol) + std::string("\n"));
    double minMixedProduct = mainModel->CheckMinimumElementMixProduct(currDispFactor);
    femUtils::WriteMessage(std::string("Minimum element Mixed Product: ") + boost::lexical_cast<std::string>(minMixedProduct) + std::string("\n"));

    if(reducedOutput){
      // Export cvPre Model ready to presolve
      std::string ncFile = "stenosis_" + std::to_string((int)currStenosisLevel) + ".coordinates";
      mainModel->WriteNodeCoordsToFile(currDispFactor,ncFile);
    }else{
      // Export cvPre Model ready to presolve
      mainModel->ExportToCvPre(currDispFactor,std::string("stenosis_") + boost::lexical_cast<std::string>(currStenosisLevel));
    }
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

// ===============
// SIMPLE MAP MODE
// ===============
int simpleMapMode(femProgramOptions* options){

  // Read All input parameters
  femInputData* data = new femInputData();
  data->ReadFromFile(options->inputFileName);

  // Read Main Model
  femModel* mainModel = new femModel();
  mainModel->ReadNodeCoordsFromFile(data->mainModelCoordsFileName,false);
  mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName,false);
  // Form Face List
  mainModel->FormElementFaceList();

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName,false);
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName,false);
  // Form Face List
  mappingModel->FormElementFaceList();

  // Read Displacement Without Rotations
  mappingModel->ReadNodeDisplacementsFromFile(data->mappingModelResultsFileName, false);

  // Export Mapping Model for Debug
  if (options->debugMode){
    mappingModel->ExportToVTKLegacy(std::string("mappingmodel.vtk"));
  }

  // Map Models
  double dispScaleFactor = 1.0;
  mainModel->MapDisplacements(mappingModel, data, dispScaleFactor);

  // Export Mapping Model for Debug
  if (options->debugMode){
    mainModel->ExportToVTKLegacy(std::string("finalmodel_simpleMap.vtk"));
  }

  // Write Application Header
  femUtils::WriteMessage(std::string("\n"));
  femUtils::WriteMessage(std::string("Completed!\n"));

  // Deallocate Pointers
  delete data;
  delete mainModel;
  delete mappingModel;

  // Done
  return 0;
}

// ========================
// TRANSLATE MODEL TO CVPRE
// ========================
int translateModelToCvPre(femProgramOptions* options){
  // Read Main Model
  femModel* mainModel = new femModel();
  bool skipFirstRow = true;

  // Read Node Coordinates and Connections
  mainModel->ReadNodeCoordsFromFile(options->inputFileName,skipFirstRow);
  mainModel->ReadElementConnectionsFromFile(options->outputFileName,skipFirstRow);

  // Form Face List
  mainModel->FormElementFaceList();

  // Export to VTK file
  mainModel->ExportToVTKLegacy(std::string("mainmodel.vtk"));

  // Orientate face nodes before exporting
  mainModel->OrientateBoundaryFaceNodes();

  // Export cvPre Model ready to presolve
  double currDispFactor = 0.0;
  mainModel->ExportToCvPre(currDispFactor,std::string("testTranslation"));

  // free heap
  delete mainModel;

  // Return
  return 0;
}

// ===============================================================
// READ MODEL AND EVALUATED VOLUME AND MIXED PRODUCT DISTRIBUTIONS
// ===============================================================
int exctractMeshQualityDistributions(int argc, char **argv){
  // Read Main Model
  femModel* mainModel = new femModel();

  // Read Node Coordinates and Connections
  mainModel->ReadNodeCoordsFromFile(std::string(argv[1]),false);
  mainModel->ReadElementConnectionsFromFile(std::string(argv[2]),false);

  // Form Box
  double limitBox[6];
  // Min x
  limitBox[0] = 58.7585;
  // Max X
  limitBox[1] = 71.8862;
  // Min Y
  limitBox[2] = 60.6504;
  // Max Y
  limitBox[3] = 70.4378;
  // Min Z
  limitBox[4] = -499.473;
  // Max Z
  limitBox[5] = -485.478;

  // Evaluate Volume and Mixed Products distributions and write it to file
  mainModel->EvalModelQualityDistributions(std::string(argv[3]),limitBox);

  // Export to VTK Legacy
  mainModel->ExportToVTKLegacy(std::string("finalmodel.vtk"));

  // Return Value
  return 0;
}

// Match The Faces in the Model
void matchModelFaces(std::string outFileName,
                     std::vector<std::string> fileList1, std::vector<femModel*> modelList1,
                     std::vector<std::string> fileList2, std::vector<femModel*> modelList2,
                     double tolerance){

  // Create Output File
  FILE* outFile;
  outFile = fopen(outFileName.c_str(),"w");

  // Loop Through all the faces
  femModel* firstModel;
  femModel* secondModel;
  for(unsigned int loopA=0;loopA<modelList1.size();loopA++){
    for(unsigned int loopB=0;loopB<modelList2.size();loopB++){
      // Store Models
      firstModel = modelList1[loopA];
      secondModel = modelList2[loopB];
      // Check Model Compatibility
      if(firstModel->isModelCompatible(secondModel,tolerance)){
        fprintf(outFile,"%s %s \n",
                femUtils::extractFileName(fileList1[loopA]).c_str(),
                femUtils::extractFileName(fileList2[loopB]).c_str());
      }
    }
  }

  // Close File
  fclose(outFile);
}

// ==================
// FIND FACE MATCHING
// ==================
int findFaceMatchList(int argc, char **argv){

  // Read the two file Names
  std::vector<std::string> fileList1;
  std::vector<std::string> fileList2;

  // Read the File Lists
  femUtils::ReadListFromFile(std::string(argv[1]),fileList1);
  femUtils::ReadListFromFile(std::string(argv[2]),fileList2);
  std::vector<femModel*> modelList1;
  std::vector<femModel*> modelList2;
  double tolerance = atof(argv[3]);

  // Fill lists
  femModel* model = nullptr;
  femUtils::WriteMessage(std::string("Reading File List 1...\n"));
  for(unsigned int loopA=0;loopA<fileList1.size();loopA++){
    model = new femModel();
   femUtils:: WriteMessage(std::string("Reading File ")+fileList1[loopA]+std::string("\n"));
    model->ReadModelNodesFromVTKFile(fileList1[loopA]);
    modelList1.push_back(model);
  }
  femUtils::WriteMessage(std::string("\n"));
  femUtils::WriteMessage(std::string("Reading File List 2...\n"));
  for(unsigned int loopA=0;loopA<fileList2.size();loopA++){
    model = new femModel();
    femUtils::WriteMessage(std::string("Reading File ")+fileList2[loopA]+std::string("\n"));
    model->ReadModelNodesFromVTKFile(fileList2[loopA]);
    modelList2.push_back(model);
  }
  // Match The Faces in the Model
  matchModelFaces(std::string("matchResult.txt"),fileList1,modelList1,fileList2,modelList2,tolerance);

  // Return Value
  return 0;
}

// ======================
// MESH VTK TRIANGULATION
// ======================
int meshVTKSkinToCVPre(int argc, char **argv){

  // Set Model File Name
  std::string VTKFileName(argv[1]);

  // Read Skin Model
  femModel* model = new femModel();
  femUtils::WriteMessage(std::string("Reading Nodes ...\n"));
  model->ReadModelNodesFromVTKFile(VTKFileName);
  femUtils::WriteMessage(std::string("Reading Elements ...\n"));
  model->ReadModelElementsFromVTKFile(VTKFileName);

  // Test: Write VTK
  model->ExportToVTKLegacy(std::string("TestExport.vtk"));

  // Convert To poly file
  std::string polyFileName("model.smesh");
  femUtils::WriteMessage(std::string("Writing SMesh File ...\n"));
  model->WriteSkinSMeshFile(polyFileName);

  // Mesh Poly file with tetgen
  femUtils::WriteMessage(std::string("Meshing with Tetgen ...\n"));
  model->MeshWithTetGen(polyFileName);

  // Export CVPRE File from node Coordinated and Element Incidences
  femUtils::WriteMessage(std::string("Exporting to CVPre ...\n"));
  translateModelToCvPre(std::string("model.1.node"),std::string("model.1.ele"),true);
}

// ==========================
// COMPUTE MODEL EXPECTATIONS
// ==========================
int computeModelExpectations(std::string modelListFileName,std::string outFile){
  // Declare Model Sequence
  femModelSequence* ms = new femModelSequence();

  // Create Model Sequence From File
  ms->ReadFromWeightedListFile(modelListFileName);

  // Compute statistics
  ms->ComputeResultStatistics();

  // Export the last Model To Vtk File
  ms->models[ms->models.size()-1]->ExportToVTKLegacy(outFile);

  // Return
  return 0;
}

// ============
// ============
// MAIN PROGRAM
// ============
// ============
int main(int argc, char **argv){

  //  Declare
  int val = 0;
  femProgramOptions* options;

  // Write Application Header
  femUtils::WriteAppHeader();

  // Get Commandline Options
  int res = options->getCommadLineOptions(argc,argv);
  if(res != 0){
    return -1;
  }

  try{
    // Normal Model Running
    switch(options->runMode){
      case rmNORMAL:
        val = runNormalMode(options);
        break;
      case rmTRANSLATETOCVPRE:
        val = translateModelToCvPre(options);
        break;
      case rmSIMPLEMAP:
        val = simpleMapMode(options);
        break;
      case rmEXTRACTMESHQUALITY:
        val = exctractMeshQualityDistributions(options);
        break;
      case rmMATCHFACELIST:
        val = findFaceMatchList(options);
        break;
      case rmMESHSKINTOCVPRE:
        val = meshVTKSkinToCVPre(options);
        break;
      case rmCOMPUTEMODELEXPECTATIONS:
        val = computeModelExpectations(options);
        break;
    }
  }catch (std::exception& ex){
    femUtils::WriteMessage(std::string(ex.what()));
    femUtils::WriteMessage(std::string("\n"));
    femUtils::WriteMessage(std::string("Program Terminated.\n"));
    return -1;
  }
  femUtils::WriteMessage(std::string("\n"));
  femUtils::WriteMessage(std::string("Program Completed.\n"));
  return val;

}
