#ifndef FEMVMSEXPLICITFLUIDSOLVER_H
#define FEMVMSEXPLICITFLUIDSOLVER_H

# include "femSolver.h"
# include "femTypes.h"
# include "femVector.h"
# include "femModel.h"

class femVMSExplicitFluidSolver: public femSolver{
  public:
    femVMSExplicitFluidSolver();
    // SOLVE PROBLEM
    virtual void solve(femModel* model);
};

#endif // FEMVMSEXPLICITFLUIDSOLVER_H
