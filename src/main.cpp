#include <stdio.h>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include "femModel.h"
#include "femInputData.h"
#include "femUtils.h"
#include "femWeightedFileName.h"
#include "femModelSequence.h"
#include "femProgramOptions.h"
#include "femSolver.h"
#include "femIncompressibleSolver.h"
#include "femSUPGExplicitFluidSolver.h"
#include "femFFDOptions.h"
#include "femPointGrid.h"
#include "femException.h"

#ifdef USE_TRILINOS

#ifdef USE_MPI
  #include "mpi.h"
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "AztecOO.h"
#include "AztecOO_config.h"
#include "Epetra_Map.h"
#include "Epetra_FEVector.h"
#include "Epetra_FECrsMatrix.h"

#endif

// ====================
// Normal Model Running
// ====================
int runNormalMode(femProgramOptions* options){
  // Vars
  double stenosisBox[6];

  // Set Debug Mode
  bool reducedOutput = options->reducedOutput;

  // Read All input parameters
  femInputData* data = new femInputData();
  data->ReadFromFile(options->inputFileName);

  // Export Main File for Debug
  if (options->debugMode){
    femUtils::ExportReferenceToVTKLegacy(data,std::string("mainreference.vtk"));
  }

  // Read Main Model
  femModel* mainModel = new femModel();
  if(!options->useVTKFile){
    // Read From TEXT Files
    mainModel->ReadNodeCoordsFromFile(data->mainModelCoordsFileName,false);
    bool numbersFromZero = false;
    mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName,false,numbersFromZero);
  }else{
    // Read From VTK
    mainModel->ReadModelNodesFromVTKFile(data->mainModelCoordsFileName);
    mainModel->ReadModelElementsFromVTKFile(data->mainModelCoordsFileName);
  }
  // Form Face List
  mainModel->FormElementFaceList();

  // Export Main File for Debug
  if (options->debugMode){
    mainModel->ExportToVTKLegacy(std::string("mainmodel.vtk"));
  }

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName,false);
  bool numbersFromZero = false;
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName,false,numbersFromZero);
  // Form Face List
  mappingModel->FormElementFaceList();

  // Read Displacement Without Rotations
  if (data->mappingDisplacementType == ipUseFile){
    mappingModel->ReadNodeDisplacementsFromFile(data->mappingModelResultsFileName, false);
  }else if(data->mappingDisplacementType == ipUseParams){
    mappingModel->ApplyParametricDisplacements(data);
  }

  // Export Mapping Model for Debug
  if (options->debugMode){
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
  if (options->debugMode){
    femUtils::ExportStenosisBoxToVTK(std::string("stenosisBox.vtk"),steNodeList);
  }  

  // Transform Mapping Model and Displacements
  femModel* newModel;
  newModel = mappingModel->TransformModel(data,stenosisBoxCenter,stenosisBox);

  // Export Mapping Model for Debug
  if (options->debugMode){
    newModel->ExportToVTKLegacy(std::string("transformedmappingmodel.vtk"));
  }

  // Map Displacement Across Models
  double dispScaleFactor = 1.0;
  mainModel->MapDisplacements(options, data, newModel, dispScaleFactor);

  // Normalize Model Displacements
  mainModel->NormalizeDisplacements(0.01*data->stenosisLength);

  // Export Mapping Model for Debug
  if (options->debugMode){
    mainModel->ExportToVTKLegacy(std::string("finalmodel.vtk"));
  }

  // Orientate face nodes before exporting
  mainModel->OrientateBoundaryFaceNodes();

  // Export cvPre Model ready to presolve
  if(!reducedOutput){
    mainModel->ExportToCvPre(0.0,std::string("reference"),options->angleLimit);
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
  if(options->debugMode){
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
    currDispFactor = mainModel->seekStenoticDisplacementFactor(data,currStenosisLevel,options->debugMode);

    // Check Volume of Undeformed Mesh
    double minVol = mainModel->CheckMinimumElementVolume(currDispFactor);
    femUtils::WriteMessage(std::string("Minimum element Volume: ") + boost::lexical_cast<std::string>(minVol) + std::string("\n"));
    double minMixedProduct = mainModel->CheckMinimumElementMixProduct(currDispFactor);
    femUtils::WriteMessage(std::string("Minimum element Mixed Product: ") + boost::lexical_cast<std::string>(minMixedProduct) + std::string("\n"));

    if(reducedOutput){
      // Export cvPre Model ready to presolve
      std::string ncFile = "stenosis_" + femUtils::intToStr((int)currStenosisLevel) + ".coordinates";
      mainModel->WriteNodeCoordsToFile(currDispFactor,ncFile);
    }else{
      // Export cvPre Model ready to presolve
      mainModel->ExportToCvPre(currDispFactor,std::string("stenosis_") + boost::lexical_cast<std::string>(currStenosisLevel),options->angleLimit);
    }
  }
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
  bool numbersFromZero = false;
  mainModel->ReadElementConnectionsFromFile(data->mainModelConnectionsFileName,false,numbersFromZero);
  // Form Face List
  mainModel->FormElementFaceList();

  // Read Mapping Model
  femModel* mappingModel = new femModel();
  mappingModel->ReadNodeCoordsFromFile(data->mappingModelCoordsFileName,false);
  mappingModel->ReadElementConnectionsFromFile(data->mappingModeConnectionslFileName,false,numbersFromZero);
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
  mainModel->MapDisplacements(options, data, mappingModel, dispScaleFactor);

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


// ===============================================================
// READ MODEL AND EVALUATED VOLUME AND MIXED PRODUCT DISTRIBUTIONS
// ===============================================================
int exctractMeshQualityDistributions(femProgramOptions* options){
  // Read Main Model
  femModel* mainModel = new femModel();

  // Read Node Coordinates and Connections
  mainModel->ReadNodeCoordsFromFile(options->inputFileName,false);
  bool numbersFromZero = false;
  mainModel->ReadElementConnectionsFromFile(options->outputFileName,false,numbersFromZero);

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
  mainModel->EvalModelQualityDistributions(std::string("graph"),limitBox);

  // Export to VTK Legacy
  mainModel->ExportToVTKLegacy(std::string("finalmodel.vtk"));

  // Return Value
  return 0;
}

// ============================
// CONVERT MODEL TO CVPRE FILES
// ============================
int translateModelToCvPre(femProgramOptions* options){
  femModel* model = new femModel();

  model->ConvertNodeAndElementsToCvPre(options->inputFileName,options->outputFileName,options->useVTKFile,false,options->angleLimit);

  delete model;

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
int findFaceMatchList(femProgramOptions* options){

  // Read the two file Names
  std::vector<std::string> fileList1;
  std::vector<std::string> fileList2;

  // Read the File Lists
  femUtils::ReadListFromFile(options->inputFileName,fileList1);
  femUtils::ReadListFromFile(options->outputFileName,fileList2);
  std::vector<femModel*> modelList1;
  std::vector<femModel*> modelList2;
  double tolerance = options->tolerance;

  // Fill lists
  femModel* model = NULL;
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
int meshVTKSkinToCVPre(femProgramOptions* options){

  // Set Model File Name
  std::string VTKFileName(options->inputFileName);

  // Read Skin Model
  femModel* model = new femModel();
  femUtils::WriteMessage(std::string("Reading Nodes ...\n"));
  model->ReadModelNodesFromVTKFile(VTKFileName);
  femUtils::WriteMessage(std::string("Reading Elements ...\n"));
  model->ReadModelElementsFromVTKFile(VTKFileName);

  // Test: Write VTK
  //model->ExportToVTKLegacy(std::string("TestExport.vtk"));

  // Convert To poly file
  std::string polyFileName("model.smesh");
  femUtils::WriteMessage(std::string("Writing SMesh File ...\n"));
  model->WriteSkinSMeshFile(polyFileName);

  // Mesh Poly file with tetgen
  femUtils::WriteMessage(std::string("Meshing with Tetgen ...\n"));
  model->MeshWithTetGen(polyFileName);

  // Export CVPRE File from node Coordinated and Element Incidences
  femUtils::WriteMessage(std::string("Exporting to CVPre ...\n"));

  // Delete Model
  delete model;

  // Create New Model
  femModel* model2 = new femModel();
  bool useVTKFile = false;
  bool skipFirstRow = true;
  model2->ConvertNodeAndElementsToCvPre(std::string("model.1.node"),std::string("model.1.ele"),useVTKFile,skipFirstRow,options->angleLimit);

  // Delete New Model
  delete model2;

  // Return OK
  return 0;
}

// ==========================
// COMPUTE MODEL EXPECTATIONS
// ==========================
int computeModelExpectations(femProgramOptions* options){
  // Declare Model Sequence
  femModelSequence* ms = new femModelSequence();

  // Create Model Sequence From File
  ms->ReadFromWeightedListFile(options->inputFileName);

  // Fix Connectivities
  // ms->FixedElementConnectivities();

  // Form Face List
  // ms->FormElementFaceList();

  // Compute wall shear stresses
  // ms->ComputeWSS();

  // Compute statistics
  ms->ComputeResultStatistics(true);

  // Export the last Model To Vtk File
  ms->models[ms->models.size()-1]->ExportToVTKLegacy(options->outputFileName);

  // Return
  return 0;
}

// =================
// COMPUTE MODEL WSS
// =================
int computeModelWSS(femProgramOptions* options){
  // Declare Model Sequence
  femModel* model = new femModel();

  // Create Model Sequence From File
  model->ReadFromVTKLegacy(options->inputFileName);

  // Write Message with total nodes and elements read
  printf("\n");
  printf("--- Total nodes found: %d\n",int(model->nodeList.size()));
  printf("--- Total elements found: %d\n",int(model->elementList.size()));
  printf("\n");

  // Fix Connectivities
  model->FixedElementConnectivities();

  // Form Face List
  model->FormElementFaceList();

  // Compute wall shear stresses
  // Units are assumed in CGS
  double viscosity = 0.04;
  model->ComputeWSS(viscosity);

  // Export to VTK
  model->ExportToVTKLegacy(options->outputFileName);

  // Return OK
  return 0;
}

// ======================
// SOLVE POISSON EQUATION
// ======================
int solvePoissonEquation(femProgramOptions* opts){

  #ifdef USE_MPI
    int *argc = NULL;
    char ***argv = NULL;
    MPI_Init(argc,argv);
  #endif

  // Create New Model
  femModel* model = new femModel();

  // Read Model From Text File
  model->ReadFromFEMTextFile(opts->inputFileName);

  // Restore Positive Volume
  model->FixedElementConnectivities();

  // Save Model For Debug
  // model->ExportToVTKLegacy("debug.vtk");

  // MPI Partitioning Information
#ifdef USE_TRILINOS
  model->CreatePartition(1);
#endif

  // CREATE NEW POISSON SOLVER
  femPoissonSolver* poisson = new femPoissonSolver();

  // CREATE OPTIONS FOR POISSON SOLVER
  femPoissonSolverOptions* slvOptions = new femPoissonSolverOptions();

  // SOLVE PROBLEM
  poisson->solve(slvOptions,model);

  // TEST : EXPORT TO VTK
  model->ExportToVTKLegacy(opts->outputFileName);

  // DELETE MODEL
  delete model;

  #ifdef USE_MPI
    MPI_Finalize();
  #endif

  // Return
  return 0;
}

#ifdef USE_TRILINOS

// ===============================================
// SOLVE STEADY-STATE ADVECTION-DIFFUSION EQUATION
// ===============================================
int solveMPISteadyStateAdvectionDiffusionEquation(femProgramOptions* options){

  // INIT MPI
  #ifdef USE_MPI
    int *argc = NULL;
    char ***argv = NULL;
    MPI_Init(argc,argv);
    Epetra_MpiComm Comm( MPI_COMM_WORLD );
  #else
    Epetra_SerialComm Comm;
  #endif
  // Print Communicator
  cout << Comm <<endl;

  // Main Node Reads Mesh and Makes Partition
  vector<femModel*> partitions;
  femModel* model;
  if(Comm.MyPID() == 0){

    // Create New Model
    model = new femModel();

    // Read Model From Text File
    model->ReadFromFEMTextFile(options->inputFileName);

    // Make Partitions    
    //partitions = model->PartitionProblem(Comm.NumProc());

    // Communicate the Models to all partitions

  }

  // Start on Current Processor
  femModel* currProcModel = model;//partitions[Comm.MyPID()];
  int NumMyElements = currProcModel->nodeList.size();

  // Construct a Map that puts same number of equations on each processor
  Epetra_Map Map(NumMyElements, 0, Comm);
  int NumGlobalElements = Map.NumGlobalElements();

  // Create a Epetra_Matrix for Finite Elements
  int rowNonZeroEst = NumGlobalElements;
  Epetra_FECrsMatrix lhs(Copy,Map,0);
  // RHS Vector
  Epetra_FEVector rhs(Map);
  // Solution Vector
  Epetra_FEVector sol(Map);

  // CREATE NEW STEADY STATE ADVECTION-DIFFUSION SOLVER
  femSteadyStateAdvectionDiffusionSolver* advDiffSolver = new femSteadyStateAdvectionDiffusionSolver();

  // CREATE OPTIONS FOR STEADY STATE ADVECTION-DIFFUSION SOLVER
  int advDiffScheme = 0;
  bool useWeakBC = 0;
  femOption* slvOptions = new femAdvectionDiffusionOptions(advDiffScheme,useWeakBC);

  // SOLVE PROBLEM
  // Assemble LHS
  advDiffSolver->assembleLHS(slvOptions,currProcModel,lhs);
  // Assemble RHS
  advDiffSolver->assembleRHS(slvOptions,currProcModel,rhs);

  // Finish up
  lhs.GlobalAssemble();
  rhs.GlobalAssemble();

  cout << lhs << endl;

  /*
  // Add Boundary Conditions
  int NumEntries = 0;
  double* Values;
  int* Indices;
  int currNode;
  double currValue;
  for(int loopA=0;loopA<model->diricheletBCNode.size();loopA++){
    currNode = model->diricheletBCNode[loopA];
    currValue = model->diricheletBCValues[loopA];
    lhs.ExtractGlobalRowCopy(currNode,0,NumEntries,Values,Indices);
    for(int loopB=0;loopB<NumEntries;loopB++){
      if(Indices[loopB] == currNode){
        Values[loopB] = 1.0;
      }else{
        Values[loopB] = 0.0;
      }
    }
    lhs.ReplaceGlobalValues(currNode,NumEntries,Values,Indices);
  }
  */

  // Apply Essential BCs to RHS
  int totDirBC = model->diricheletBCNode.size();
  int GIDs[totDirBC];
  double GIDVal[totDirBC];
  for(int loopA=0;loopA<totDirBC;loopA++){
    GIDs[loopA] = model->diricheletBCNode[loopA];
    GIDVal[loopA] = model->diricheletBCValues[loopA][0];
  }
  rhs.ReplaceGlobalValues(totDirBC,GIDs,GIDVal);

  // Create Linear Problem
  Epetra_LinearProblem problem(&lhs, &sol, &rhs);
   // Create AztecOO instance
  AztecOO solver(problem);

  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.Iterate(10, 1.0E-8);

  solver.PrintLinearSystem("systemPrint.out");

  cout << "Solver performed " << solver.NumIters() << " iterations." << endl
       << "Norm of true residual = " << solver.TrueResidual() << endl;

  double theNorm = 0.0;
  (void) rhs.Norm2 (&theNorm);
  cout << "Norm of sol: " << theNorm << endl;

  cout << sol << endl;

  #ifdef USE_MPI
    MPI_Finalize() ;
  #endif

return 0;

}

#endif

// ===============================================
// SOLVE STEADY-STATE ADVECTION-DIFFUSION EQUATION
// ===============================================
int solveSteadyStateAdvectionDiffusionEquation(femProgramOptions* options){
  // Create New Model
  femModel* model = new femModel();

  // Read Model From Text File
  model->ReadFromFEMTextFile(options->inputFileName);

  // CREATE NEW STEADY STATE ADVECTION-DIFFUSION SOLVER
  femSteadyStateAdvectionDiffusionSolver* advDiffSolver = new femSteadyStateAdvectionDiffusionSolver();

  // CREATE OPTIONS FOR POISSON SOLVER

  int advDiffScheme = 0;
  bool useWeakBC = true;
  femOption* slvOptions = new femAdvectionDiffusionOptions(advDiffScheme,useWeakBC);

  // SOLVE PROBLEM
  advDiffSolver->solve(slvOptions,model);

  // Export Model to Check
  model->ExportToVTKLegacy(string("out.vtk"));

  // Return
  return 0;
}

// =======================
// SOLVE INCOMPRESSIBLE NS
// =======================
int solveIncompressibleNS(femProgramOptions* options){

  // Create New Model
  femModel* model = new femModel();

  // Read Model From Text File
  model->ReadFromFEMTextFile(options->inputFileName);

  // CREATE NEW STEADY STATE ADVECTION-DIFFUSION SOLVER
  femIncompressibleSolver* incNS = new femIncompressibleSolver();

  // SOLVE PROBLEM
  incNS->solve(model);

  // Return
  return 0;
}

// ========================
// TEST ELEMENT FORMULATION
// ========================
int testElementFormulation(femProgramOptions* options){
  // Create New Model
  femModel* model = new femModel();

  // Read Node Coord
  bool skipFirstRow = false;
  model->ReadNodeCoordsFromFile(options->nodeFileName,skipFirstRow);

  // Read Element Connections
  bool numbersFromZero = true;
  model->ReadElementConnectionsFromFile(options->connectionFileName,skipFirstRow,numbersFromZero);

  // Restore Positive Volume
  model->FixedElementConnectivities();

  // CREATE NEW POISSON SOLVER
  femTestSolver* test = new femTestSolver();

  // CREATE OPTIONS FOR POISSON SOLVER
  femTestSolverOptions* slvOptions = new femTestSolverOptions();

  // SOLVE PROBLEM
  test->solve(slvOptions,model);

  // Return
  return 0;
}

// =========================
// RUN FREE-FORM DEFORMATION
// =========================
int runFFD(femProgramOptions* options){

  // Read FFD Input File
  femFFDOptions* ffdOpts = new femFFDOptions();
  ffdOpts->readFromFile(options->inputFileName);

  // Read VTK Model
  femModel* model = new femModel();

  // Create Model Sequence From File
  model->ReadFromVTKLegacy(ffdOpts->inputFileName);

  // Create the vector of FFD PointGrids
  std::vector<femPointGrid*> ffdGrids;
  for(int loopA=0;loopA<ffdOpts->ffdData.size();loopA++){

    // Create a new FFD Grid with the data from file
    femPointGrid* ffdGrid = new femPointGrid(ffdOpts->ffdData[loopA].minPoint,
                                             ffdOpts->ffdData[loopA].totPoints,
                                             ffdOpts->ffdData[loopA].gridAxis_S,
                                             ffdOpts->ffdData[loopA].gridAxis_T,
                                             ffdOpts->ffdData[loopA].gridAxis_U,
                                             ffdOpts->ffdData[loopA].dispNodes,
                                             ffdOpts->ffdData[loopA].dispVals);

    // Add to the collection of FFD Transformations
    ffdGrids.push_back(ffdGrid);
  }

  // Apply the deformations on the model
  string outGrid("bernsteinGrid_");
  for(int loopA=0;loopA<ffdGrids.size();loopA++){
    ffdGrids[loopA]->morphModel(model);
    ffdGrids[loopA]->ExportToVTKLegacy(outGrid + to_string(loopA+1) + string(".vtk"));
  }

  // Export Model to VTK
  model->ExportToVTKLegacy(options->outputFileName);

  // Return OK
  return 0;
}

// =============================
// SOLVE EXPLICIT NAVIER-STOCKES 
// =============================
int solveExplicitNS(femProgramOptions* options){

  // Create New Model
  femModel* model = new femModel();

  // Read Model From Text File
  model->ReadFromFEMTextFile(options->inputFileName);

  // CREATE NEW STEADY STATE ADVECTION-DIFFUSION SOLVER
  // femVMSExplicitFluidSolver* exNS = new femVMSExplicitFluidSolver();
  femSUPGExplicitFluidSolver* exNS = new femSUPGExplicitFluidSolver();

  // SOLVE PROBLEM
  exNS->solve(model);

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
  femProgramOptions* options = new femProgramOptions();

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
      case rmCOMPUTEMODELWSS:
        val = computeModelWSS(options);
        break;
      case rmSOLVEPOISSON:
        val = solvePoissonEquation(options);
        break;
      case rmSOLVESTEADYSTATEADVECTIONDIFFUSION:
        val = solveSteadyStateAdvectionDiffusionEquation(options);
        //val = solveMPISteadyStateAdvectionDiffusionEquation(options);
        break;
      case rmSOLVEINCOMPRESSIBLENS:
        val = solveIncompressibleNS(options);
        break;
      case rmTESTELEMENTS:
        val = testElementFormulation(options);
        break;
      case rmFFD:
        val = runFFD(options);
        break;
      case rmEXPLICITNS:
        val = solveExplicitNS(options);
        break;        
    }
  }catch(femException ex){
    femUtils::WriteMessage(std::string(ex.what()));
    femUtils::WriteMessage(std::string("\n"));
    femUtils::WriteMessage(std::string("Program Terminated.\n"));
    femUtils::WriteMessage(std::string("\n"));
    return -1;
  }
  femUtils::WriteMessage(std::string("\n"));
  femUtils::WriteMessage(std::string("Program Completed.\n"));
  femUtils::WriteMessage(std::string("\n"));
  delete options;
  return val;

}
