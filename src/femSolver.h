#ifndef FEMSOLVER_H
#define FEMSOLVER_H

# include "femModel.h"
# include "femOption.h"
# include "femMatrix.h"
# include "femVector.h"

#ifdef USE_TRILINOS
# include "femTrilinosMatrix.h"
# include "femTrilinosVector.h"
# include "Epetra_FEVector.h"
# include "Epetra_SerialDenseVector.h"
# include "Epetra_FECrsMatrix.h"
#endif


// GENERIC SOLVER
class femSolver{
  public:
    femSolver();
    // SOLVE PROBLEM
    virtual void       solve(femOption* options, femModel* model);
#ifdef USE_ARMADILLO
    femDoubleVec       solveLinearSystem(femDenseMatrix* poissonMat,femDenseVector* poissonVec);
#endif
#ifdef USE_CSPARSE
    femDoubleVec       solveLinearSystem(femSparseMatrix* lhs,femDenseVector* rhs);
#endif
#ifdef USE_TRILINOS
    femTrilinosVector* solveLinearSystem(femTrilinosMatrix* lhs,femTrilinosVector* rhs);
#endif
};

// STEADY STATE ADVECTION DIFFUSION SOLVER
class femSteadyStateAdvectionDiffusionSolver: public femSolver{
  public:
    femSteadyStateAdvectionDiffusionSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
    // Additional Routines
#ifdef USE_TRILINOS
    void assembleLHS(femOption* options, femModel* model,Epetra_FECrsMatrix &lhs);
    void assembleRHS(femOption* options, femModel* model,Epetra_FEVector &rhs);
#endif

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
