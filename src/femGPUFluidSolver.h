#ifndef FEMGPUFLUIDSOLVER_H
#define FEMGPUFLUIDSOLVER_H

# include "femSolver.h"
# include "femTypes.h"
# include "femVector.h"
# include "femModel.h"

class femGPUFluidSolver: public femSolver{
  public:
    femGPUFluidSolver();
    // SOLVE PROBLEM
    virtual void solve(femModel* model);
};

#endif // FEMGPUFLUIDSOLVER_H
