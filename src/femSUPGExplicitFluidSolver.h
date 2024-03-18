#ifndef FEMSUPGEXPLICITFLUIDSOLVER_H
#define FEMSUPGEXPLICITFLUIDSOLVER_H

# include "femSolver.h"
# include "femTypes.h"
# include "femVector.h"
# include "femModel.h"

class femSUPGExplicitFluidSolver: public femSolver{
  public:
    femSUPGExplicitFluidSolver();
    // SOLVE PROBLEM
    virtual void solve(femModel* model);
};

#endif // FEMSUPGEXPLICITFLUIDSOLVER_H
