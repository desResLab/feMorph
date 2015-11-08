# include "femVector.h"
# include "femModel.h"

// ==========
// FEM VECTOR
// ==========

// CONSTRUCTOR
femVector::femVector(){
  //throw femException("femVector constructor not implemented.\n");
}

// GET VECTOR SIZE
int femVector::getSize(){
  throw femException("femVector getSize not implemented.\n");
}

// GET VECTOR COMPONENT
double femVector::getComponent(int id){
  throw femException("femVector getComponent not implemented.\n");
}

// SET VECTOR COMPONENT
double femVector::setComponent(int id, double entry){
  throw femException("femVector setComponent not implemented.\n");
}

// ASSEMBLE IN DENSE VECTOR
void femVector::assemble(femDoubleVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}

// ASSEMBLE IN DENSE VECTOR
void femVector::blockAssemble(femDoubleBlockVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}

// APPLY DIRICHELET VECTORS
void femVector::applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues){
  throw femException("Not Implemented.\n");
}

// WRITE VECTOR TO FILE
void femVector::writeToFile(string fileName){
  throw femException("Not Implemented.\n");
}

// ================
// FEM DENSE VECTOR
// ================

// CONSTRUCTOR
femDenseVector::femDenseVector(int total){
  values.resize(total);
  for(int loopA=0;loopA<total;loopA++){
    values[loopA] = 0.0;
  }
}

// GET VECTOR COMPONENT
double femDenseVector::getComponent(int id){
  return values[id];
}

// SET VECTOR COMPONENT
double femDenseVector::setComponent(int id, double entry){
  values[id] = entry;
}

// ASSEMBLE IN DENSE VECTOR
void femDenseVector::blockAssemble(femDoubleBlockVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}

// ASSEMBLE IN DENSE VECTOR
void femDenseVector::assemble(femDoubleVec vec,femIntVec indices){
  int currIndex = 0;
  for(size_t loopA=0;loopA<indices.size();loopA++){
    currIndex = indices[loopA];
    values[currIndex] += vec[loopA];
  }
}

// APPLY DIRICHELET VECTORS
void femDenseVector::applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues){
  for(size_t loopA=0;loopA<diricheletBCNode.size();loopA++){
    values[diricheletBCNode[loopA]] = diricheletBCValues[loopA];
  }
}

// WRITE VECTOR TO FILE
void femDenseVector::writeToFile(string fileName){
  //Create File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  for(size_t loopA=0;loopA<values.size();loopA++){
    fprintf(f,"%e \n",values[loopA]);
  }
  // Close File
  fclose(f);
}

void femVector::GlobalAssemble(){
  throw femException("femVector::GlobalAssemble Not Implemented.\n");
}

void femDenseVector::GlobalAssemble(){
  throw femException("femDenseVector::GlobalAssemble Not Implemented.\n");
}
