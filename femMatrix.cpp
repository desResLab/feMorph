# include "femMatrix.h"
# include "femException.h"

femMatrix::femMatrix(){
  throw femException("Not Implemented.\n");
}

// CONSTRUCTOR FOR DENSE MATRIX
femDenseMatrix::femDenseMatrix(femModel* model){
  int totDof = model->nodeList.size();
  values.resize(totDof);
  for(int loopA=0;loopA<totDof;loopA++){
    values[loopA].resize(totDof);
  }
  for(int loopA=0;loopA<totDof;loopA++){
    for(int loopB=0;loopB<totDof;loopB++){
      values[loopA][loopB] = 0.0;
    }
  }
}

// CONSTRUCTOR FOR SPARSE MATRIX
femSparseMatrix::femSparseMatrix(femModel* model){
  throw femException("Not Implemented.\n");
}

// FEM MATRIX ASSEMBLE: NOT IMPLEMENTED
void femMatrix::assemble(femDoubleMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}

// ASSEMBLE DENSE MATRIX
void femDenseMatrix::assemble(femDoubleMat elMat,femIntVec connections){
  int rowIndex = 0;
  int colIndex = 0;
  int count = elMat.size();
  for(int loopA=0;loopA<count;loopA++){
    rowIndex = connections[loopA];
    for(int loopB=0;loopB<count;loopB++){
      colIndex = connections[loopB];
      values[rowIndex][colIndex] += elMat[loopA][loopB];
    }
  }
}

// ASSEMBLE SPARSE MATRIX
void femSparseMatrix::assemble(femDoubleMat elMat,femIntVec connections){
  int rowIndex = 0;
  int colIndex = 0;
  int counter = 0;
  bool found = false;
  int count = elMat.size();
  for(int loopA=0;loopA<count;loopA++){
    rowIndex = connections[loopA];
    for(int loopB=0;loopB<count;loopB++){
      colIndex = connections[loopB];
      // Find the Start of the column
      counter = diagPtr[colIndex];
      found = false;
      while((!found)&&(counter<diagPtr[colIndex+1]-1)){
        found = (rowIndex == rowPtr[counter]);
        // Update Counter
        if(!found){
          counter++;
        }
      }
      if(!found){
        throw femException("ERROR: femSparseMatrix::assemble NOT FOUND.\n");
      }else{
        values[counter] += elMat[loopA][loopB];
      }
    }
  }
}

// MAIN CLASS: NOT IMPLEMENTED
void femMatrix::applyDirichelet(femIntVec dofs){
  throw femException("Not Implemented.\n");
}

// APPLY DIRICHELET CONDITION TO DENSE MATRIX
// ALL ELEMENTS ARE ZERO AND ONE IN THE DIAGONAL
void femDenseMatrix::applyDirichelet(femIntVec dofs){
  int size = dofs.size();
  int currIndex = 0;
  for(int loopA=0;loopA<size;loopA++){
    currIndex = dofs[loopA];
    // Operate on row
    for(int loopB=0;loopB<(int)values.size();loopB++){
      if(loopB == currIndex){
        values[currIndex][loopB] = 1.0;
      }else{
        values[currIndex][loopB] = 0.0;
      }
    }
    // Operate on column
    for(int loopB=0;loopB<(int)values[0].size();loopB++){
      if(loopB == currIndex){
        values[loopB][currIndex] = 1.0;
      }else{
        values[loopB][currIndex] = 0.0;
      }
    }
  }
}

// APPLY DIRICHELET CONDITIONS TO SPARSE MATRIX
void femSparseMatrix::applyDirichelet(femIntVec dofs){
  int size = dofs.size();
  int currIndex = 0;
  for(int loopA=0;loopA<size;loopA++){
    currIndex = dofs[loopA];
    // Transverse the matrix and set to zero the elements
    for(int loopB=0;loopB<totCols;loopB++){
      // If Current Column, delete all value except diagonal
      if(loopB == currIndex){
        values[diagPtr[loopB]] = 1.0;
        for(int loopC=diagPtr[loopB]+1;loopC<diagPtr[loopB+1];loopC++){
          values[loopC] = 0.0;
        }
      }else{
        for(int loopC=diagPtr[loopB];loopC<diagPtr[loopB+1];loopC++){
          if(rowPtr[loopC] == currIndex){
            values[loopC] = 0.0;
          }
        }
      }
    }
  }
}




