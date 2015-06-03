#ifndef FEMSOLVER_H
#define FEMSOLVER_H

# include "femModel.h"
# include "femOption.h"
# include "femMatrix.h"
# include "femVector.h"
# include "Epetra_FEVector.h"
# include "Epetra_SerialDenseVector.h"
# include "Epetra_FECrsMatrix.h"


// GENERIC SOLVER
class femSolver{
  public:
    femSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
    femDoubleVec solveLinearSystem(femSparseMatrix* lhs,femVector* rhs);
    femDoubleVec solveLinearSystem(femDenseMatrix* poissonMat,femVector* poissonVec);
};

// STEADY STATE ADVECTION DIFFUSION SOLVER
class femSteadyStateAdvectionDiffusionSolver: public femSolver{
  public:
    femSteadyStateAdvectionDiffusionSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
    // Additional Routines
    void assembleLHS(femOption* options, femModel* model,Epetra_FECrsMatrix &lhs);
    void assembleRHS(femOption* options, femModel* model,Epetra_FEVector &rhs);

};

// POISSON SOLVER
class femPoissonSolver: public femSolver{
  public:
    femPoissonSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
};

// TEST ELEMENTS
class femTestSolver: public femSolver{
  public:
    femTestSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
};


#endif // FEMSOLVER_H
