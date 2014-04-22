#ifndef FEMSOLVER_H
#define FEMSOLVER_H

#include "femModel.h"

// GENERIC SOLVER
class femSolver{
  public:
    femSolver();
    // SOLVE PROBLEM
    //virtual void solve(femModel* model);
};

// ADVECTION DIFFUSION SOLVER
class femAdvectionDiffusionSolver: public femSolver{
  public:
    femAdvectionDiffusionSolver();
    // SOLVE PROBLEM
    //virtual void solve(femAdvectionDiffusionOptions* options, femModel* model);
};


#endif // FEMSOLVER_H
