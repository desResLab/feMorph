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
};

// DERIVED DENSE MATRIX CLASS
class femDenseMatrix: public femMatrix{
  public:
    // Matrix Values
    femDoubleMat values;

    // CONSTRUCTOR
    femDenseMatrix(femModel* model);

    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
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

    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
};


#endif // FEMMATRIX_H
