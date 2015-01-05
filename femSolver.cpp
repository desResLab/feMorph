# include <armadillo>

# include "femSolver.h"
# include "femMatrix.h"
# include "femVector.h"
# include "femException.h"

// MAIN SOLVER CLASS
femSolver::femSolver(){
}

// ADVECTION-DIFFUSION SOLVER
femAdvectionDiffusionSolver::femAdvectionDiffusionSolver(){
}

// POISSON SOLVER
femPoissonSolver::femPoissonSolver(){
}

// GENERIC SOLVER: NOT IMPLEMENTED
void femSolver::solve(femOption* options, femModel* model){
  throw femException("Not Implemented.\n");
}

// SOLVE ADVECTION-DIFFUSION PROBLEM
void femAdvectionDiffusionSolver::solve(femOption* options, femModel* model){
/*
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

  }*/
}

// ======================================
// SOLVE LINEAR SYSTEM WITH DIRECT METHOD
// ======================================
femDoubleVec solveLinearSystem(femDenseMatrix* poissonMat,femVector* poissonVec){
  // Initialize Armadillo Matrix
  arma::mat A(poissonMat->totRows,poissonMat->totCols);
  for(int loopA=0;loopA<poissonMat->totRows;loopA++){
    for(int loopB=0;loopB<poissonMat->totCols;loopB++){
      A(loopA,loopB) = poissonMat->values[loopA][loopB];
    }
  }
  // Fill RHS Vector
  arma::vec b(poissonVec->values.size());
  for(int loopA=0;loopA<poissonVec->values.size();loopA++){
    b(loopA) = poissonVec->values[loopA];
  }
  // Solve Linear Set of equations
  arma::vec x = solve(A, b);
  // Return
  femDoubleVec result;
  for(int loopA=0;loopA<x.size();loopA++){
    result.push_back(x[loopA]);
  }
}

// ======================
// SOLVE POISSON EQUATION
// ======================
void femPoissonSolver::solve(femOption* options, femModel* model){

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irFirstOrder);

  // Init Sparse Matrix and Dense Vector
  femMatrix* poissonMat;
  poissonMat = new femDenseMatrix(model);
  femVector* poissonVec = new femVector((int)model->nodeList.size());

  // Local Element Matrix
  femDoubleMat elMat;
  femDoubleVec elSourceVec;
  femDoubleVec elBCVec;
  femIntVec currIdxList;

  // ASSEMBLE MATRIX FROM ALL ELEMENTS
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Gauss Points Loop
    // Assemble Poisson Matrix
    model->elementList[loopElement]->formPoissonMatrix(model->nodeList,rule,elMat);
    // Assemble Sparse Matrix
    poissonMat->assemble(elMat,model->elementList[loopElement]->elementConnections);
  }

  // ASSEMBLE SOURCE TERM
  int currEl = 0;
  double currValue = 0.0;
  for(size_t loopSource=0;loopSource<model->sourceElement.size();loopSource++){
    // Get Current Element
    currEl = model->sourceElement[loopSource];
    // Eval Source Vector
    model->elementList[currEl]->formPoissonSource(model->nodeList,rule,model->sourceValues[loopSource],elSourceVec);
    // Assemble Source Vector
    poissonVec->assemble(elSourceVec,model->elementList[currEl]->elementConnections);
  }

  // ASSEMBLE NEUMANN BOUNDARY CONDITIONS
  for(size_t loopBC=0;loopBC<model->neumannBCElement.size();loopBC++){
    // Get Current Element
    currEl = model->neumannBCElement[loopBC][0];
    currValue = model->neumannBCValues[loopBC];
    // Assemble Local Indexes
    currIdxList.clear();
    for(size_t loopA=1;loopA<model->neumannBCElement.size();loopA++){
      currIdxList.push_back(model->neumannBCElement[loopBC][loopA]);
    }
    // Eval Boundary Flux Vector
    model->elementList[currEl]->formPoissonBCFlux(model->nodeList,currIdxList,currValue,rule,elBCVec);
    // Assemble in RHS
    poissonVec->assemble(elBCVec,currIdxList);
  }

  // APPLY DIRICHELET BOUNDARY CONDITIONS
  // Sparse Matrix
  poissonMat->applyDirichelet(model->diricheletBCNode);
  // RHS Vector
  poissonVec->applyDirichelet(model->diricheletBCNode,model->diricheletBCValues);

  // SOLVE LINEAR SYSTEM OF EQUATIONS
  femDoubleVec solution;
  solution = solveLinearSystem((femDenseMatrix*)poissonMat,poissonVec);
}

