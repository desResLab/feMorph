#include "femOption.h"

femOption::femOption(){
}

femAdvectionDiffusionOptions::femAdvectionDiffusionOptions(int scheme, string fileName){
  advDiffScheme = scheme;
  outputFileName = fileName;
}

femPoissonSolverOptions::femPoissonSolverOptions(){
}

femTestSolverOptions::femTestSolverOptions(){

}

