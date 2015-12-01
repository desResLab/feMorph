#ifndef FEMSOLVER_H
#define FEMSOLVER_H

# include "femModel.h"
# include "femOption.h"
# include "femMatrix.h"
# include "femVector.h"
# include "femTime.h"

#ifdef USE_TRILINOS
# include "trilinos/femTrilinosMatrix.h"
# include "trilinos/femTrilinosVector.h"
# include "Epetra_FEVector.h"
# include "Epetra_SerialDenseVector.h"
# include "Epetra_FECrsMatrix.h"
#endif


// GENERIC SOLVER
class femSolver{
  public:
    femSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
#ifdef USE_ARMADILLO
    femVector* solveLinearSystem(femDenseMatrix* poissonMat,femDenseVector* poissonVec);
#endif
#ifdef USE_CSPARSE
    femVector* solveLinearSystem(femSparseMatrix* lhs,femDenseVector* rhs);
#endif
#ifdef USE_TRILINOS
    // SOLVE LINEAR SYSTEM
    femVector* solveLinearSystem(int totalNodes, femTrilinosMatrix* lhs, femTrilinosVector* rhs, int nodeDOFs);
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
