#include "femOption.h"

femOption::femOption(){
}

femAdvectionDiffusionOptions::femAdvectionDiffusionOptions(int scheme, string fileName, bool locUseWeakBC){
  advDiffScheme = scheme;
  outputFileName = fileName;
  useWeakBC = locUseWeakBC;
}

femPoissonSolverOptions::femPoissonSolverOptions(){
}

femTestSolverOptions::femTestSolverOptions(){

}

