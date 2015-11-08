#ifndef FEMMATRIX_H
#define FEMMATRIX_H

# include "femModel.h"
# include "femTypes.h"

// MAIN MATRIX CLASS
class femMatrix{
  public:
    // CONSTRUCTOR
    femMatrix();

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void blockAssemble(femDoubleBlockMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
    virtual void clearRowAndColumn(int dof);
    virtual void completeFill();
};

// DERIVED DENSE MATRIX CLASS
class femDenseMatrix: public femMatrix{
  public:
    // Store Rows and Columns
    int totRows;
    int totCols;

    // Matrix Values
    femDoubleMat values;

    // CONSTRUCTOR
    femDenseMatrix(femModel* model);

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void blockAssemble(femDoubleBlockMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
    virtual void clearRowAndColumn(int dof);
    virtual void completeFill();
};


// DERIVED SPARSE MATRIX CLASS
class femSparseMatrix: public femMatrix{
  public:
    // Store Columns Only
    int totCols;

    // Matrix Values
    femIntVec diagPtr; // Pointer to Columns
    femIntVec rowPtr; // Row indices
    femDoubleVec values; // Numerical Values

    // CONSTRUCTOR
    femSparseMatrix(femModel* model);

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void blockAssemble(femDoubleBlockMat elMat,femIntVec connections);

    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
    virtual void clearRowAndColumn(int dof);
    virtual void completeFill();
};


#endif // FEMMATRIX_H
