#include "femOption.h"

femOption::femOption(){
}

femAdvectionDiffusionOptions::femAdvectionDiffusionOptions(int scheme, int velType, int sourceType, string fileName){
  advDiffScheme = scheme;
  advDiffVelType = velType;
  advDiffSourceType = sourceType;
  outputFileName = fileName;
}

femPoissonSolverOptions::femPoissonSolverOptions(){
}

femTestSolverOptions::femTestSolverOptions(){

}

