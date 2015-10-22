#include "femTrilinosVector.h"

// CONSTRUCTOR
femTrilinosVector::femTrilinosVector(int total,int nodeDOFs){
  values.resize(total);
  for(int loopA=0;loopA<total;loopA++){
    values[loopA].resize(nodeDOFs);
    for(int loopB=0;loopB<total;loopB++){
      values[loopA][loopB] = 0.0;
    }
  }
}

// GETTER AND SETTER
int femTrilinosVector::getSize(){
  throw femException("Not Implemented.\n");
}
// ASSEMBLE IN DENSE COLUMN FORMAT
void femTrilinosVector::assemble(femDoubleVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}
void femTrilinosVector::assembleDOF(femDoubleDOFVec vec,femIntVec indices){
  throw femException("Not Implemented.\n");
}
void femTrilinosVector::applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues){
  throw femException("Not Implemented.\n");
}
void femTrilinosVector::writeToFile(string fileName){
  throw femException("Not Implemented.\n");
}
