# include "femMatrix.h"
# include "femException.h"

femMatrix::femMatrix(){
}

// CONSTRUCTOR FOR DENSE MATRIX
femDenseMatrix::femDenseMatrix(femModel* model){
  int totDof = model->nodeList.size();
  totRows = totDof;
  totCols = totDof;
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

// =============================
// CONSTRUCTOR FOR SPARSE MATRIX
// =============================
femSparseMatrix::femSparseMatrix(femModel* model){

    // Get Total Number of Equations
  int totDof = model->nodeList.size();

  // Assign Rows and Columns
  totRows = totDof;
  totCols = totDof;

  // Form Column Pointer
  // Initialize
  diagPtr.resize(totDof+1);
  for(int loopA=0;loopA<(totDof+1);loopA++){
    diagPtr[loopA] = 0;
  }
  // Fill node Counter
  int nodeCounter[totDof];
  for(int loopA=0;loopA<totDof;loopA++){
    nodeCounter[loopA] = 1;
  }
  int currSize = 0.0;
  int currNode = 0;
  for(size_t loopA=0;loopA<model->elementList.size();loopA++){
    currSize = model->elementList[loopA]->elementConnections.size();
    for(size_t loopB=0;loopB<currSize;loopB++){
      currNode = model->elementList[loopA]->elementConnections[loopB];
      nodeCounter[currNode] += currSize - 1;
    }
  }

  // Form Column Pointer
  int totNoZero = 0;
  for(int loopA=1;loopA<totDof;loopA++){
    diagPtr[loopA] = diagPtr[loopA-1] + nodeCounter[loopA];
    totNoZero += nodeCounter[loopA];
  }
  diagPtr[totDof] = totNoZero;

  // Initialize Values
  values.resize(totNoZero);
  for(int loopA=0;loopA<totNoZero;loopA++){
    values[loopA] = 0.0;
  }

  // Form Pointer to Row Elements
  femDoubleMat tempRowPtrMat;
  tempRowPtrMat.resize(totDof);
  // Fill with Diagonal Elements
  for(int loopA=0;loopA<totDof;loopA++){
    tempRowPtrMat[loopA].push_back(loopA);
  }
  // Fill with extra diagonal elements
  int currNode1 = 0;
  int currNode2 = 0;
  for(size_t loopA=0;loopA<model->elementList.size();loopA++){
    for(size_t loopB=0;loopB<model->elementList[loopA]->elementConnections.size();loopB++){
      currNode1 = model->elementList[loopA]->elementConnections[loopB];
      for(size_t loopC=0;loopC<model->elementList[loopA]->elementConnections.size();loopC++){
        currNode2 = model->elementList[loopA]->elementConnections[loopC];
        if(currNode1 != currNode2){
          tempRowPtrMat[currNode1].push_back(currNode2);
        }
      }
    }
  }
  // Sort All Entries
  for(int loopA=0;loopA<totDof;loopA++){
    std::sort(tempRowPtrMat[loopA].begin(), tempRowPtrMat[loopA].end());
  }

  // Copy to Row Pointer
  rowPtr.resize(totNoZero);
  int count = 0;
  for(int loopA=0;loopA<totDof;loopA++){
    for(int loopB=0;loopB<tempRowPtrMat[loopA].size();loopB++){
      rowPtr[count] = tempRowPtrMat[loopA][loopB];
      count++;
    }
  }

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
      while((!found)&&(counter<diagPtr[colIndex+1])){
        found = (rowIndex == rowPtr[counter]);
        // Update Counter
        if(!found){
          counter++;
        }
      }
      if(!found){
        printf("RowIndex %d, ColIndex %d\n",rowIndex,colIndex);
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
    for(int loopB=0;loopB<(int)values[0].size();loopB++){
      if(loopB == currIndex){
        values[currIndex][loopB] = 1.0;
      }else{
        values[currIndex][loopB] = 0.0;
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

// WRITE MATRIX TO FILE
void femMatrix::writeToFile(string fileName){
  throw femException("Not Implemented.\n");
}

void femDenseMatrix::writeToFile(string fileName){
  //Create File
  FILE* f;
  f = fopen(fileName.c_str(),"w");
  for(int loopA=0;loopA<totRows;loopA++){
    for(int loopB=0;loopB<totCols;loopB++){
      fprintf(f,"%e ",values[loopA][loopB]);
    }
    fprintf(f,"\n");
  }
  // Close File
  fclose(f);
}

void femSparseMatrix::writeToFile(string fileName){
  throw femException("Not Implemented.\n");
}

double femMatrix::getRowSum(int loopA){
  throw femException("Not Implemented.\n");
}

double femDenseMatrix::getRowSum(int row){
  double sum = 0.0;
  for(int loopA=0;loopA<totCols;loopA++){
    sum += values[row][loopA];
  }
  return sum;
}

double femSparseMatrix::getRowSum(int loopA){
  throw femException("Not Implemented.\n");
}



