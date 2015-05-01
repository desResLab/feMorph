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
  int advDiffVelType;
  int advDiffSourceType;
  string outputFileName;
  // Constructor
  femAdvectionDiffusionOptions(int scheme, int velType, int sourceType, string fileName);
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
