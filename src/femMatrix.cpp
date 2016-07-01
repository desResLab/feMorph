# include <algorithm>

# include "femMatrix.h"

femMatrix::femMatrix(){
}

// ============================
// CONSTRUCTOR FOR DENSE MATRIX
// ============================
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

    double getModelNodalTopology(femIntVec& diagPtr,femIntVec& rowPtr);
  // Get Total Number of Equations
  int totDof = model->nodeList.size();

  // Assign Rows and Columns
  totCols = totDof;

  // Form Pointer to Row Elements
  femIntMat tempRowPtrMat;
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
          if(find(tempRowPtrMat[currNode1].begin(), tempRowPtrMat[currNode1].end(), currNode2) == tempRowPtrMat[currNode1].end()){
            tempRowPtrMat[currNode1].push_back(currNode2);
          }
        }
      }
    }
  }

  // Sort All Entries
  for(int loopA=0;loopA<totDof;loopA++){
    std::sort(tempRowPtrMat[loopA].begin(), tempRowPtrMat[loopA].end());
  }
  // Print Index Array
  //printf("INDEX ARRAY\n");
  //for(int loopA=0;loopA<totDof;loopA++){
  //  for(int loopB=0;loopB<tempRowPtrMat[loopA].size();loopB++){
  //    printf("%d ",tempRowPtrMat[loopA][loopB]);
  //  }
  //  printf("\n");
  //}

  // Initialize Column Pointer
  diagPtr.resize(totDof+1);
  for(int loopA=0;loopA<(totDof+1);loopA++){
    diagPtr[loopA] = 0;
  }

  // Form Column Pointer
  int totNoZero = tempRowPtrMat[0].size();
  for(int loopA=1;loopA<totDof;loopA++){
    diagPtr[loopA] = diagPtr[loopA-1] + tempRowPtrMat[loopA-1].size();
    totNoZero += tempRowPtrMat[loopA].size();
  }
  diagPtr[totDof] = totNoZero;

  // Initialize Values
  values.resize(totNoZero);
  for(int loopA=0;loopA<totNoZero;loopA++){
    values[loopA] = 0.0;
  }

  // Copy to Row Pointer
  rowPtr.resize(totNoZero);
  int count = 0;
  for(int loopA=0;loopA<totDof;loopA++){
    for(size_t loopB=0;loopB<tempRowPtrMat[loopA].size();loopB++){
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
  int count = connections.size();
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
      for(int loopC=diagPtr[loopB];loopC<diagPtr[loopB+1];loopC++){
        if(rowPtr[loopC] == currIndex){
          if(currIndex == loopB){
            values[loopC] = 1.0;
          }else{
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

// PLOT SPARSE MATRIX TO FILE
void femSparseMatrix::writeToFile(std::string fileName){
  //Create File
  FILE* f;
  f = fopen(fileName.c_str(),"w");

  int totDof = totCols;

  // Plot Diagonal Pointer
  fprintf(f,"# Column Pointer\n");
  for(size_t loopA=0;loopA<diagPtr.size();loopA++){
    fprintf(f,"%d \n",diagPtr[loopA]);
  }

  // Plot Row Index Pointer
  fprintf(f,"# Row Index Pointer\n");
  for(int loopA=0;loopA<totDof;loopA++){
    for(int loopB=diagPtr[loopA];loopB<diagPtr[loopA+1];loopB++){
      fprintf(f,"%d ",rowPtr[loopB]);
    }
    fprintf(f,"\n");
  }

  // Plot Matrix Values
  fprintf(f,"# Matrix Values\n");
  for(size_t loopA=0;loopA<values.size();loopA++){
    fprintf(f,"%f \n",values[loopA]);
  }

  // Close File
  fclose(f);
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

// SET ROW AND COLUMN TO ZERO IN DENSE MATRIX
void femMatrix::clearRowAndColumn(int dof){
  throw femException("Not Implemented.\n");
}

// SET ROW AND COLUMN TO ZERO IN DENSE MATRIX
void femDenseMatrix::clearRowAndColumn(int dof){
  for(int loopA=0;loopA<totRows;loopA++){
    values[loopA][dof] = 0.0;
  }
  for(int loopA=0;loopA<totCols;loopA++){
    values[dof][loopA] = 0.0;
  }
}

// SET ROW AND COLUMN TO ZERO IN SPARSE MATRIX
void femSparseMatrix::clearRowAndColumn(int dof){
  // Set to Zero in Column dof
  for(int loopA=diagPtr[dof];loopA<diagPtr[dof+1];loopA++){
    values[loopA] = 0.0;
  }
  // Set to Zero with dof Row
  for(size_t loopA=0;loopA<rowPtr.size();loopA++){
    if(rowPtr[loopA] == dof){
      values[loopA] = 0.0;
    }
  }
}

// Block Assemble only implemented for Trilinos Matrices
void femMatrix::blockAssemble(femDoubleBlockMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}
void femDenseMatrix::blockAssemble(femDoubleBlockMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}
void femSparseMatrix::blockAssemble(femDoubleBlockMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}

// Complete Fill only implemented for Trilinos Matrices
void femMatrix::completeFill(){
  throw femException("femMatrix::completeFill not Implemented.\n");
}
void femDenseMatrix::completeFill(){
  throw femException(" femDenseMatrix::completeFill not Implemented.\n");
}
void femSparseMatrix::completeFill(){
  throw femException("femSparseMatrix::completeFill not Implemented.\n");
}
