#ifndef FEMOPTION_H
#define FEMOPTION_H

# include <string>

# include "femTypes.h"

using namespace std;

// GENERAL FEM OPTION
class femOption{
  public:
    femOption();
};

// ADVECTION DIFFUSION SOLVER OPTIONS
class femAdvectionDiffusionOptions: public femOption{
public:
  // Data Members
  int advDiffScheme;
  bool useWeakBC;
  // Constructor
  femAdvectionDiffusionOptions(int scheme, bool useWeakBC);
};

// TRANSIENT ADVECTION DIFFUSION SOLVER OPTIONS
class femIncompressibleSolverOption: public femOption{
public:
  // Data Members
  int nodeDofs;
  // Time Integration Parameters
  double timeStep;
  int totalTimeSteps;
  double alphaM;
  double alphaF;
  double gamma;
  // Stages of Integration 0-NS, other AdvectionDiffusion
  femIntVec stages;
  // Use Prescribed Velocities or Solve them from NS
  bool usePrescribedVelocity;
  int prescribedVelType;
  // Constructor
  femIncompressibleSolverOption();
};

// POISSON SOLVER OPTIONS
class femPoissonSolverOptions: public femOption{
  public:
    femPoissonSolverOptions();
};

// TEST SOLVER OPTIONS
class femTestSolverOptions: public femOption{
  public:
    femTestSolverOptions();
};

#endif // FEMOPTION_H
