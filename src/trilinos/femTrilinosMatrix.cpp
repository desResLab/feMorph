# include "femTrilinosMatrix.h"

// CONSTRUCTOR
femTrilinosMatrix::femTrilinosMatrix(femModel* model,int nodeDOFs){
  // Create Communicator
  #ifdef USE_MPI
    EpetraComm = new Epetra_MpiComm(MPI_COMM_WORLD);
  #else
    EpetraComm = new Epetra_SerialComm();
  #endif

  // Set Node DOFs
  this->nodeDOFs = nodeDOFs;

  // Need to create the Block matrix first
  map  = new Epetra_BlockMap(-1,
                             model->totNodesInProc,
                             model->localToGlobalNodes,
                             nodeDOFs,
                             0,
                             *EpetraComm);
  // Create Graph
  graph = new Epetra_CrsGraph(Copy,*map,0);

  // Get Model Topology
  femIntVec diagPtr;
  femIntVec rowPtr;
  model->getModelNodalTopology(diagPtr,rowPtr);

  // Fill Graph Indices
  int NumIndices = 0;
  for(int loopA=0;loopA<model->totNodesInProc;loopA++){
    NumIndices = diagPtr[loopA+1] - diagPtr[loopA];
    graph->InsertGlobalIndices(loopA,NumIndices, &rowPtr[diagPtr[loopA]]);
  }
  graph->FillComplete();

  // Create the Epetra_FEVector
  values = new Epetra_FEVbrMatrix(Copy,*graph);
}

femTrilinosMatrix::~femTrilinosMatrix(){
  delete map;
  delete graph;
  delete values;
}

// VIRTUAL FUNCTIONS
void  femTrilinosMatrix::assemble(femDoubleMat nodeMat,femIntVec elConnections){
  if(nodeDOFs > 1){
    throw femException("ERROR: Calling Assemble with more than one DOF per node.\n");
  }
  // Convert Vector to int*
  int currNodeRow = 0;
  int currNodeCol = 0;
  // Loop on the Block Entries
  // Loop on Row
  for(size_t loopA=0;loopA<elConnections.size();loopA++){
    currNodeRow = elConnections[loopA];
    // Start Summing Values in epetra block FE Matrix
    values->BeginSumIntoGlobalValues(currNodeRow,(int)elConnections.size(),&elConnections[0]);
    // Lop on Column
    for(size_t loopB=0;loopB<elConnections.size();loopB++){
      currNodeCol = elConnections[loopA];
      // Assign Entries to this block
      values->SubmitBlockEntry(&nodeMat[loopA][loopB],nodeDOFs,nodeDOFs,nodeDOFs);
    }
    // Finish submitting block row entries
    values->EndSubmitEntries();
  }
  // Global Assemble
  values->GlobalAssemble();
}
void  femTrilinosMatrix::blockAssemble(femDoubleBlockMat nodeMat, femIntVec elConnections){
  // Convert Vector to int*
  int currNodeRow = 0;
  int currNodeCol = 0;
  // Loop on the Block Entries
  // Loop on Row
  for(size_t loopA=0;loopA<elConnections.size();loopA++){
    currNodeRow = elConnections[loopA];
    // Start Summing Values in epetra block FE Matrix
    values->BeginSumIntoGlobalValues(currNodeRow,(int)elConnections.size(),&elConnections[0]);
    // Lop on Column
    for(size_t loopB=0;loopB<elConnections.size();loopB++){
      currNodeCol = elConnections[loopA];
      // Assign Entries to this block
      values->SubmitBlockEntry(&nodeMat[0][loopA][loopB],nodeDOFs,nodeDOFs,nodeDOFs);
    }
    // Finish submitting block row entries
    values->EndSubmitEntries();
  }
  // Global Assemble
  values->GlobalAssemble();
}

// SET ALL THE EXTRADIAGONALS TO 0.0 AND DIAGONALS TO 1.0
void  femTrilinosMatrix::applyBlockDirichelet(femIntVec gNodesIdx,int dof){

  // Convert Vector to int*
  int currNodeRow = 0;
  int currNodeCol = 0;
  int rowDOFs = 0;
  int NumColumnBlockEntries = 0;
  int* ColBlockIndices = NULL;
  int* ColDims = NULL;
  double newValues[nodeDOFs*nodeDOFs];
  int test = 0;

  // Get Max Number of Block Columns Entries
  int MaxNumBlockEntries = values->GlobalMaxNumBlockEntries();

  // Allocate Column Indexes and Dimensions
  ColBlockIndices = new int[MaxNumBlockEntries];
  ColDims         = new int[MaxNumBlockEntries];

  // Loop on the Block Entries
  // Loop on Row
  for(size_t loopA=0;loopA<gNodesIdx.size();loopA++){

    // Set Current Row
    currNodeRow = gNodesIdx[loopA];

    // Start Extracting Global Rows - Get the indices of the available column blocks
    test = values->BeginExtractGlobalBlockRowCopy(currNodeRow,MaxNumBlockEntries,rowDOFs,NumColumnBlockEntries,ColBlockIndices,ColDims);

    // Start Inserting Values in epetra block FE Matrix
    test = values->BeginReplaceGlobalValues(currNodeRow,NumColumnBlockEntries,ColBlockIndices);

    // Loop on Column
    for(size_t loopB=0;loopB<NumColumnBlockEntries;loopB++){
      // Get current global node index
      currNodeCol = ColBlockIndices[loopB];

      if(currNodeRow != currNodeCol){
        // Extradiagonal should be zero
        for(int loopC=0;loopC<nodeDOFs;loopC++){
          for(int loopD=0;loopD<nodeDOFs;loopD++){
            newValues[loopC*nodeDOFs + loopD] = 0.0;
          }
        }
      }else{
        // Diagonal Block
        for(int loopC=0;loopC<nodeDOFs;loopC++){
          for(int loopD=0;loopD<nodeDOFs;loopD++){
            if(loopC == loopD){
              newValues[loopC*nodeDOFs + loopD] = 1.0;
            }else{
              newValues[loopC*nodeDOFs + loopD] = 0.0;
            }
          }
        }
      }
      // Assign Entries to this block
      values->SubmitBlockEntry(newValues,nodeDOFs,nodeDOFs,nodeDOFs);
    }
    // Finish submitting block row entries
    values->EndSubmitEntries();
  }
  // Free
  delete [] ColBlockIndices;
  delete [] ColDims;
}

// ======================================================
// APPLY DIRICHELET CONDITIONS ON MATRIX: NOT IMPLEMENTED
// ======================================================
void  femTrilinosMatrix::applyDirichelet(femIntVec dofs){
  throw femException("Not Implemented.\n");
}

// WRITE MATRIX TO FILE
void   femTrilinosMatrix::writeToFile(string fileName){
  // Open File
  ofstream myfile(fileName.c_str());
  // Print Vector
  values->Print (myfile);
  // Close Output File
  myfile.close();
}
double femTrilinosMatrix::getRowSum(int loopA){
  throw femException("Not Implemented.\n");
}
void   femTrilinosMatrix::clearRowAndColumn(int dof){
  throw femException("Not Implemented.\n");
}

// COMPLETE FILLING OF MATRIX
void femTrilinosMatrix::completeFill(){
  values->FillComplete();
}
