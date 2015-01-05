# include "femVector.h"
# include "femModel.h"

// CONSTRUCTOR
femVector::femVector(int total){
  values.resize(total);
  for(int loopA=0;loopA<total;loopA++){
    values[loopA] = 0.0;
  }
}

// ASSEMBLE IN DENSE VECTOR
void femVector::assemble(femDoubleVec vec,femIntVec indices){
  int currIndex = 0;
  for(size_t loopA=0;loopA<indices.size();loopA++){
    currIndex = indices[loopA];
    values[currIndex] += vec[loopA];
  }
}

// APPLY DIRICHELET VECTORS
void femVector::applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues){
  for(int loopA=0;loopA<diricheletBCNode.size();loopA++){
    values[diricheletBCNode[loopA]] = diricheletBCValues[loopA];
  }
}
