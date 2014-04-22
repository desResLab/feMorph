#include "femSolver.h"

femSolver::femSolver(){
}

femAdvectionDiffusionSolver::femAdvectionDiffusionSolver(){
}


// SOLVE ADVECTION-DIFFUSION PROBLEM
/*void femAdvectionDiffusionSolver::solve(femAdvectionDiffusionOptions* options, femModel* model){

  // Initialize Solution
  initSolution(model,scalarVector);

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irFirstOrder);

  // Time Loop
  double currentTime = 0.0;
  for(int loopTime=0;loopTime<options.totalTimeSteps;loopTime++){

    // Advance in time
    if(loopTime>0){
      currentTime += options.timeStep;
    }

    // Get Local Mass and Stiffness Matrices from each element
    for(int loopElement=0;loopElement<model->elementList.size();loopElement++){
      // Get Stabilitzation Parameter at Gauss Points
      elementList[loopElement]->assembleStab(rule,tauSUPG);
      // Get Mass
      elementList[loopElement]->assembleMass(nodeVelocities,model->nodeList,tauSUPG,rule,massMat);
      // Get Stiffness
      elementList[loopElement]->assembleStiffness(nodeVelocities,model->nodeList,tauSUPG,rule,diffusivity,stiffnessMat);
      // Get RHS
      elementList[loopElement]->assembleRHS();
      // Assemble in Global Sparse Matrix and RHS Vector
      assembleInSparseMatrixAndRHS(massMat,stiffnessMat,RHS);
    }

    // Solve Linear System

  }
}*/

