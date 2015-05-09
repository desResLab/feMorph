#include "femOption.h"

femOption::femOption(){
}

femAdvectionDiffusionOptions::femAdvectionDiffusionOptions(int scheme, bool locUseWeakBC){
  advDiffScheme = scheme;
  useWeakBC = locUseWeakBC;
}

femPoissonSolverOptions::femPoissonSolverOptions(){
}

femTestSolverOptions::femTestSolverOptions(){

}

