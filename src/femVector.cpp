# include "femVector.h"
# include "femModel.h"

// ==========
// FEM VECTOR
// ==========

// CONSTRUCTOR
femVector::femVector(){
  throw femException("Not Implemented.\n");
}

// ASSEMBLE IN DENSE VECTOR
void femVector::assemble(femDoubleVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}

// ASSEMBLE IN DENSE VECTOR
void femVector::assembleDOF(femDoubleDOFVec vec,femIntVec indices){
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

// ASSEMBLE IN DENSE VECTOR
void femDenseVector::assembleDOF(femDoubleDOFVec vec,femIntVec indices){
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
