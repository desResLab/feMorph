#ifndef FEMSOLVER_H
#define FEMSOLVER_H

# include "femModel.h"
# include "femOption.h"

// GENERIC SOLVER
class femSolver{
  public:
    femSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
};

// ADVECTION DIFFUSION SOLVER
class femAdvectionDiffusionSolver: public femSolver{
  public:
    femAdvectionDiffusionSolver();
    // SOLVE PROBLEM
    virtual void solve(femOption* options, femModel* model);
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
