#ifndef FEMOPTION_H
#define FEMOPTION_H

// GENERAL FEM OPTION
class femOption{
  public:
    femOption();
};

// ADVECTION DIFFUSION SOLVER OPTIONS
class femAdvectionDiffusionOptions: public femOption{
public:
  femAdvectionDiffusionOptions();
};

// POISSON SOLVER OPTIONS
class femPoissonSolverOptions: public femOption{
  public:
    femPoissonSolverOptions();
};

#endif // FEMOPTION_H
