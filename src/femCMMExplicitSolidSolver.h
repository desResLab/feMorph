#ifndef FEMCMMEXPLICITSOLIDSOLVER_H
#define FEMCMMEXPLICITSOLIDSOLVER_H

# include "femSolver.h"
# include "femTypes.h"
# include "femVector.h"
# include "femModel.h"

class femCMMExplicitSolidSolver: public femSolver{
  public:
    femCMMExplicitSolidSolver();
    // SOLVE PROBLEM
    virtual void solve(femModel* model);
};

#endif // FEMCMMEXPLICITSOLIDSOLVER_H
