# include <string>
# include <armadillo>

# include "femSolver.h"
# include "femMatrix.h"
# include "femVector.h"
# include "femException.h"

using namespace std;

// MAIN SOLVER CLASS
femSolver::femSolver(){
}

// ADVECTION-DIFFUSION SOLVER
femAdvectionDiffusionSolver::femAdvectionDiffusionSolver(){
}

// POISSON SOLVER
femPoissonSolver::femPoissonSolver(){
}

// TEST SOLVER
femTestSolver::femTestSolver(){

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
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

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
  printf("Solving...\n");
  // Solve Linear Set of equations
  arma::vec x = solve(A, b);
  // Return
  femDoubleVec result;
  for(int loopA=0;loopA<x.size();loopA++){
    result.push_back(x[loopA]);
  }
  // Return Result
  return result;
}

// ======================
// SOLVE POISSON EQUATION
// ======================
void femPoissonSolver::solve(femOption* options, femModel* model){

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Init Sparse Matrix and Dense Vector
  femMatrix* poissonMat;
  poissonMat = new femDenseMatrix(model);
  femVector* poissonVec = new femVector((int)model->nodeList.size());

  // Local Element Matrix
  femDoubleMat elMat;
  femDoubleVec elSourceVec;
  femDoubleVec elBCVec;
  femIntVec currIdxList;

  // Print One local Matrix
  //double sum = 0.0;
  //printf("Local Matrix\n");
  //model->elementList[0]->formPoissonMatrix(model->nodeList,rule,elMat);
  //Create File
  //FILE* f;
  //f = fopen("pLocalMat.txt","w");
  //for(int loopA=0;loopA<elMat.size();loopA++){
  //  for(int loopB=0;loopB<elMat.size();loopB++){
  //    fprintf(f,"%e ",elMat[loopA][loopB]);
  //  }
  //  fprintf(f,"\n");
  //}
  // Close File
  //fclose(f);

  // ASSEMBLE MATRIX FROM ALL ELEMENTS
  printf("Assembling Matrix...\n");
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Gauss Points Loop
    // Assemble Poisson Matrix
    model->elementList[loopElement]->formPoissonMatrix(model->nodeList,rule,(femDoubleVec)model->elDiffusivity[loopElement],elMat);
    // Assemble Sparse Matrix
    poissonMat->assemble(elMat,model->elementList[loopElement]->elementConnections);
  }

  // ASSEMBLE SOURCE TERM
  printf("Assembling Source...\n");
  int currEl = 0;
  double currDiff = 0.0;
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
  printf("Assembling Neumann...\n");
  int currNode = 0;
  elBCVec.resize(model->nodeList.size());
  for(int loopA=0;loopA<model->nodeList.size();loopA++){
    elBCVec[loopA] = 0.0;
  }
  for(size_t loopBC=0;loopBC<model->neumannBCElement.size();loopBC++){
    // Get Current Element
    currEl = model->neumannBCElement[loopBC];
    // Get Value
    currValue = model->neumannBCValues[loopBC];
    // Get Global Face Nodes
    for(int loopA=0;loopA<model->neumannBCFaceNodes[loopBC].size();loopA++){
      currNode = model->neumannBCFaceNodes[loopBC][loopA];
      currDiff = model->elDiffusivity[currEl][0];
      elBCVec[currNode] += (currDiff * currValue)/(double)model->neumannBCFaceNodes[loopBC].size();
    }
  }
  // Assemble in RHS
  for(int loopA=0;loopA<poissonVec->getSize();loopA++){
    poissonVec->values[loopA] += elBCVec[loopA];
  }


  // APPLY DIRICHELET BOUNDARY CONDITIONS
  printf("Assembling Dirichelet...\n");
  // Sparse Matrix
  poissonMat->applyDirichelet(model->diricheletBCNode);
  // RHS Vector
  poissonVec->applyDirichelet(model->diricheletBCNode,model->diricheletBCValues);

  // PRINT MATRIX AND VECTOR TO FILE
  //poissonMat->writeToFile(string("pMatrix.txt"));
  //poissonVec->writeToFile(string("pVector.txt"));

  // Check the sum of terms in the matrix
  //sum = 0.0;
  //for(int loopA=0;loopA<poissonMat->totRows;loopA++){
  //  sum = poissonMat->getRowSum(loopA);
  //  printf("SUM: %e\n",sum);
  //}

  // SOLVE LINEAR SYSTEM OF EQUATIONS
  femDoubleVec solution;
  solution = solveLinearSystem((femDenseMatrix*)poissonMat,poissonVec);

  // ADD SOLUTION TO MODEL RESULTS
  femResult* res = new femResult();
  res->label = string("PoissonResult");
  res->type = frNode;
  // Assign to values
  for(size_t loopA=0;loopA<solution.size();loopA++){
    //printf("Solution %d %f\n",loopA,solution[loopA]);
    res->values.push_back(solution[loopA]);
  }
  model->resultList.push_back(res);

  // Free Memory
  //delete poissonMat;
  //delete poissonVec;
}

// ========================
// TEST ELEMENT FORMULATION
// ========================
void femTestSolver::solve(femOption* options, femModel* model){

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Create a Unit Nodal Quantity
  femDoubleVec unitNodalFunction;
  unitNodalFunction.resize(model->nodeList.size());
  double xCoord,yCoord,zCoord;
  for(int loopA=0;loopA<model->nodeList.size();loopA++){
    // Get Coordinates for the current Node
    xCoord = model->nodeList[loopA]->coords[0];
    yCoord = model->nodeList[loopA]->coords[1];
    zCoord = model->nodeList[loopA]->coords[2];
    unitNodalFunction[loopA] = 1.0;
  }

  // Compute the total Volume using Gaussian Integration
  double volume = 0.0;
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Assemble Poisson Matrix
    volume += model->elementList[loopElement]->integrateNodalVector(model->nodeList,rule,unitNodalFunction);
  }
  printf("Total Volume: %f\n",volume);

  // Create a Nodal Quantity
  femDoubleVec nodalFunction;
  nodalFunction.resize(model->nodeList.size());
  for(int loopA=0;loopA<model->nodeList.size();loopA++){
    // Get Coordinates for the current Node
    xCoord = model->nodeList[loopA]->coords[0];
    yCoord = model->nodeList[loopA]->coords[1];
    zCoord = model->nodeList[loopA]->coords[2];
    nodalFunction[loopA] = 2.0*xCoord -1.0*yCoord + zCoord;
  }

  // Compute the Gradient of this Quantity at various locatios
  int currElement = 0;
  int numNodes = model->elementList[currElement]->numberOfNodes;
  // Get Gauss Points and Weights
  femDoubleMat intCoords = rule->getCoords(numNodes);

  // Eval Current Shape Derivatives Matrix at the first Gauss Point
  femDoubleMat shapeDeriv;
  model->elementList[currElement]->evalGlobalShapeFunctionDerivative(model->nodeList,intCoords[0][0],intCoords[0][1],intCoords[0][2],shapeDeriv);

  // Compute the Gradient
  femDoubleVec grad;
  int currNode;
  grad.resize(3);
  for(int loopA=0;loopA<kDims;loopA++){
    grad[loopA] =  0.0;
    for(int loopB=0;loopB<numNodes;loopB++){
      currNode = model->elementList[currElement]->elementConnections[loopB];
      grad[loopA] += shapeDeriv[loopB][loopA] * nodalFunction[currNode];
    }
  }

  // Print Gradient
  printf("Gradient X: %f\n",grad[0]);
  printf("Gradient Y: %f\n",grad[1]);
  printf("Gradient Z: %f\n",grad[2]);
}

