#ifndef FEMOPTION_H
#define FEMOPTION_H

# include <string>

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
