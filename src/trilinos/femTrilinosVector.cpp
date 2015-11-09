# include "femTrilinosVector.h"

// CONSTRUCTOR
femTrilinosVector::femTrilinosVector(int totNodesInProc,int* localToGlobal,int nodeDOFs){
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
                             totNodesInProc,
                             localToGlobal,
                             nodeDOFs,
                             0,
                             *EpetraComm);

  // Create the Epetra_FEVector
  //values = new Epetra_FEVector(*map,totNodesInProc);
  values = new Epetra_FEVector(*map,nodeDOFs);
}

// CONSTRUCTOR
femTrilinosVector::femTrilinosVector(Epetra_BlockMap map,int totComponents,int nodeDOFs){
  // Set Node DOFs
  this->nodeDOFs = nodeDOFs;

  // Create the Epetra_FEVector
  //values = new Epetra_FEVector(map,totComponents);
    values = new Epetra_FEVector(map,nodeDOFs);
}

// GETTER AND SETTER
int femTrilinosVector::getSize(){
  return values->GlobalLength();
}

// GET VECTOR COMPONENT
double femTrilinosVector::getComponent(int id){
  int currBlock = id / nodeDOFs;
  int currIndex = id % nodeDOFs;
  double** pointer;
  values->ExtractView(&pointer);
  // Get Data
  return pointer[currIndex][currBlock];
}

// SET VECTOR COMPONENT
void femTrilinosVector::setComponent(int id, double entry){
  int currBlock = id / nodeDOFs;
  int currIndex = id % nodeDOFs;
  double** pointer;
  values->ExtractView(&pointer);
  // Get Data
  pointer[currIndex][currBlock] = entry;
}

// ADD COMPONENT TO A VECTOR
/*void femTrilinosVector::addComponent(int id, double entry){
  int currBlock = id / nodeDOFs;
  int currIndex = id % nodeDOFs;
  double** pointer;
  double value = 0.0;
  values->ExtractView(&pointer);
  // Get Data
  value = pointer[currIndex][currBlock];
  pointer[currIndex][currBlock] = value + entry;
}*/

void femTrilinosVector::addComponent(int id, double entry){
  int GIDs[1];
  GIDs[0] = id;
  double val[1];
  val[0] = entry;
  values->SumIntoGlobalValues(1,GIDs,val);
  values->GlobalAssemble();
}

// ASSEMBLE VECTOR INDICES WITH 1 DOF
void femTrilinosVector::assemble(femDoubleVec vec,femIntVec indices){
  if(nodeDOFs > 1){
    throw femException("ERROR: Calling Assemble with more than one DOF per node.\n");
  }
  /*int idx[1];
  double nums[1];
  for(int loopA=0;loopA<(int)indices.size();loopA++){
    idx[0] = indices[loopA];
    nums[0] = vec[loopA];
    values->SumIntoGlobalValues(1,idx,nums);
  }*/
  values->SumIntoGlobalValues(indices.size(),&indices[0],&vec[0]);
  values->GlobalAssemble();
}

// ASSEMBLE BLOCK VECTOR
void femTrilinosVector::blockAssemble(femDoubleBlockVec vec,femIntVec indices){
  int currIdx = 0.0;
  double currVal = 0.0;
  // Loop through the Nodes
  for(size_t loopA=0;loopA<indices.size();loopA++){
    currIdx = indices[loopA];
    // Loop through the entries for every node
    for(size_t loopB=0;loopB<vec[loopA].size();loopB++){
      currVal = vec[loopA][loopB];
      values->SumIntoGlobalValue(currIdx,loopB,currVal);
    }
  }
}

// APPLY DIRICHELET BOUNDARY CONDITIONS
void femTrilinosVector::applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues){
  int currIdx = 0.0;
  double currVal = 0.0;
  // Loop through the Nodes
  for(size_t loopA=0;loopA<diricheletBCNode.size();loopA++){
    currIdx = diricheletBCNode[loopA];
    currVal = diricheletBCValues[loopA];
    // Loop through the entries for every node
    for(size_t loopB=0;loopB<nodeDOFs;loopB++){
      values->ReplaceGlobalValue(currIdx,loopB,currVal);
    }
  }
}

// WRITE VECTOR TO FILE
void femTrilinosVector::writeToFile(string fileName){
  // Open File
  ofstream myfile(fileName.c_str());
  // Print Vector
  values->Print (myfile);
  // Close Output File
  myfile.close();
}

// COMPLETE ASSEMBLY
void femTrilinosVector::GlobalAssemble(){
  values->GlobalAssemble();
}

