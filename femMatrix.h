#ifndef FEMMATRIX_H
#define FEMMATRIX_H

# include "femModel.h"
# include "femTypes.h"

// MAIN MATRIX CLASS
class femMatrix{
  public:
    int totRows;
    int totCols;

    // CONSTRUCTOR
    femMatrix();

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
};

// DERIVED DENSE MATRIX CLASS
class femDenseMatrix: public femMatrix{
  public:
    // Matrix Values
    femDoubleMat values;

    // CONSTRUCTOR
    femDenseMatrix(femModel* model);

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);

};


// DERIVED SPARSE MATRIX CLASS
class femSparseMatrix: public femMatrix{
  public:
    // Matrix Values
    femIntVec diagPtr; // Pointer to Columns
    femIntVec rowPtr; // Row indices
    femDoubleVec values; // Numerical Values

    // CONSTRUCTOR
    femSparseMatrix(femModel* model);

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
};


#endif // FEMMATRIX_H