#ifndef FEMINCOMPRESSIBLESOLVER_H
#define FEMINCOMPRESSIBLESOLVER_H

# include "femSolver.h"

class femIncompressibleSolver: public femSolver{
  public:
    femIncompressibleSolver();
    // SOLVE PROBLEM
    virtual void solve(femModel* model);
};

#endif // FEMINCOMPRESSIBLESOLVER_H
