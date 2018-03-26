# include <string>

#ifdef USE_ARMADILLO
  # include <armadillo>
#endif

# include "femSolver.h"
# include "femMatrix.h"
# include "femVector.h"
# include "femException.h"

#ifdef USE_CSPARSE
extern "C" {
#include "csparse.h"
}
#endif

#ifdef USE_TRILINOS
  # include "Epetra_LinearProblem.h"
  # include "AztecOO.h"
#endif

using namespace std;

#ifdef USE_ARMADILLO
// =========================================================
// SOLVE LINEAR SYSTEM WITH DIRECT METHOD WITH ARMADILLO LIB
// =========================================================
femDoubleVec femSolver::solveLinearSystem(femDenseMatrix* poissonMat,femDenseVector* poissonVec){
  // Initialize Armadillo Matrix
  arma::mat A(poissonMat->totRows,poissonMat->totCols);
  for(int loopA=0;loopA<poissonMat->totRows;loopA++){
    for(int loopB=0;loopB<poissonMat->totCols;loopB++){
      A(loopA,loopB) = poissonMat->values[loopA][loopB];
    }
  }
  // Fill RHS Vector
  arma::vec b(poissonVec->values.size());
  for(size_t loopA=0;loopA<poissonVec->values.size();loopA++){
    b(loopA) = poissonVec->values[loopA];
  }
  printf("Solving...\n");
  // Solve Linear Set of equations
  arma::vec x = arma::solve(A, b);
  // Return
  femDoubleVec result;
  for(int loopA=0;loopA<x.size();loopA++){
    result.push_back(x[loopA]);
  }
  // Return Result
  return result;
}
#endif

#ifdef USE_CSPARSE
// =============================================
// SOLVE SPARSE LINEAR SYSTEM WITH DIRECT METHOD
// =============================================
femVector* femSolver::solveLinearSystem(femSparseMatrix* lhs,femDenseVector* rhs){
  // Set Tolerance
  double tol = 1.0e-12;
  // Set Sparse Ordering AMD
  int sparseOrdering = 1;
  // Init Matrix
  cs* A;
  // Fill Matrix in Compressed Column Mode
  int totDofs = lhs->totCols;
  int totNonZeros = lhs->rowPtr.size();
  A->nzmax = totNonZeros;
  A->m = lhs->totCols;
  A->n = lhs->totCols;

  // Assign Matrix Structure
  A->p = new int[totDofs+1];
  double* b = new double[totDofs];
  A->i = new int[totNonZeros];
  A->x = new double[totNonZeros];
  for(int loopA=0;loopA<totDofs+1;loopA++){
    A->p[loopA] = lhs->diagPtr[loopA];
  }
  for(int loopA=0;loopA<totNonZeros;loopA++){
    A->i[loopA] = lhs->rowPtr[loopA];
    A->x[loopA] = lhs->values[loopA];
  }
  A->nz = -1;
  // Copy rhs in PCorr
  for(int loopA=0;loopA<totDofs;loopA++){
    b[loopA] = rhs->values[loopA];
  }
  // Print RHS
  // printf("RHS Before Solution\n");
  // for(int loopA=0;loopA<totalNodes;loopA++){
  //   printf("%d %e\n",loopA,pCorr[loopA]);
  // }

  // Solve system

  int ok = cs_lusol(A,b,sparseOrdering,tol);
  if (ok == 0){
    std::string errorMsg("Error: Cannot Solve Linear System\n");
    throw femException(errorMsg.c_str());
  }
  // Print Solution
  // Copy rhs in PCorr
  // printf("Pressure Corrections\n");
  // for(int loopA=0;loopA<totalNodes;loopA++){
  //   printf("%d %e\n",loopA,pCorr[loopA]);
  // }

  // Copy Solution Back
  femVector* result = new femDenseVector(totDofs);
  for(int loopA=0;loopA<totDofs;loopA++){
    result->setComponent(loopA,b[loopA]);
  }
  // Return Result
  return result;
}
#endif

#ifdef USE_TRILINOS
// ======================================
// SOLVE LINEAR SYSTEM WITH DIRECT METHOD
// ======================================
femVector* femSolver::solveLinearSystem(int totalNodes, femTrilinosMatrix* lhs,femTrilinosVector* rhs,int nodeDOFs){

  femTrilinosVector* sol = new femTrilinosVector(lhs->values->RowMap(),totalNodes,nodeDOFs);

  Epetra_LinearProblem problem(&*lhs->values,&*sol->values,&*rhs->values);

  // Create AztecOO instance
  AztecOO solver(problem);

  solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  //solver.SetAztecOption(AZ_subdomain_solve, AZ_lu);
  //solver.SetAztecOption(AZ_precond, AZ_ls);

  solver.Iterate(1000, 1.0E-5);

  const double* status = solver.GetAztecStatus();

  cout << "Solver performed " << solver.NumIters() << " iterations." << endl
  << "Norm of true residual = " << solver.TrueResidual() << endl;

  // Return Solution
  return sol;
}
#endif

// MAIN SOLVER CLASS
femSolver::femSolver(){
}

// ADVECTION-DIFFUSION SOLVER
femSteadyStateAdvectionDiffusionSolver::femSteadyStateAdvectionDiffusionSolver(){
}

// POISSON SOLVER
femPoissonSolver::femPoissonSolver(){
}

// TEST SOLVER
femTestSolver::femTestSolver(){

}

// GENERIC SOLVER: NOT IMPLEMENTED
void femSolver::solve(femOption* options, femModel* model){
  throw femException("Solve Not Implemented for femSolver.\n");
}

// =================================
// SOLVE ADVECTION-DIFFUSION PROBLEM
// =================================
void femSteadyStateAdvectionDiffusionSolver::solve(femOption* options, femModel* model){
    printf("Solving Advection-Diffusion...\n");

    // One DOF per Node in Advection Diffusion
    int nodeDOFs = 1;

    // Get Integration Rule
    femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

    // Init Sparse Matrix and Dense Vector
    femMatrix* advDiffMat;
    femVector* advDiffVec;
#ifdef USE_ARMADILLO
    advDiffMat = new femDenseMatrix(model);
    advDiffVec = new femDenseVector((int)model->nodeList.size());
#endif
#ifdef USE_CSPARSE
    advDiffMat = new femSparseMatrix(model);
    advDiffVec = new femDenseVector((int)model->nodeList.size());
#endif
#ifdef USE_TRILINOS
    advDiffMat = new femTrilinosMatrix(model,nodeDOFs);
    // One Nodal Degrees of freedom in AdvDiff
    advDiffVec = new femTrilinosVector(model->totNodesInProc,model->localToGlobalNodes,nodeDOFs);
#endif

    // Local Element Matrix
    femDoubleMat elMat;
    femDoubleVec elRhs;

    // ASSEMBLE LHS MATRIX
    printf("Assembling LHS Matrix...\n");
    for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
      // Assemble Advection-Diffusion Matrix
      model->elementList[loopElement]->formAdvDiffLHS(model->nodeList,
                                                      rule,
                                                      (femDoubleVec)model->elDiffusivity[loopElement],
                                                      (femDoubleVec)model->elVelocity[loopElement],
                                                      ((femAdvectionDiffusionOptions*)options)->advDiffScheme,
                                                      elMat);
      // Assemble Sparse Matrix
      advDiffMat->assemble(elMat,model->elementList[loopElement]->elementConnections);
    }

    // ASSEMBLE RHS TERM
    printf("Assembling RHS Vector...\n");
    int currEl = 0;
    double currDiff = 0.0;
    double currValue = 0.0;
    for(size_t loopSource=0;loopSource<model->sourceElement.size();loopSource++){
      // Get Current Element
      currEl = model->sourceElement[loopSource];

      model->sourceValues[loopSource],
      (femDoubleVec)model->elDiffusivity[loopSource],
      (femDoubleVec)model->elVelocity[loopSource],
      ((femAdvectionDiffusionOptions*)options)->advDiffScheme,

      // Eval Source Vector
      model->elementList[currEl]->formAdvDiffRHS(model->nodeList,
                                                 rule,
                                                 model->sourceValues[loopSource],(femDoubleVec)model->elDiffusivity[loopSource],(femDoubleVec)model->elVelocity[loopSource],
                                                 ((femAdvectionDiffusionOptions*)options)->advDiffScheme,
                                                 elRhs);
      // Assemble Source Vector
      advDiffVec->assemble(elRhs,model->elementList[currEl]->elementConnections);
    }

    // APPLY DIRICHELET/ESSENTIAL BOUNDARY CONDITIONS

    if(((femAdvectionDiffusionOptions*)options)->useWeakBC){
      // USE WEAK BOUNDARY CONDITIONS
      printf("Assembling Weak Boundary Conditions...\n");

      // REMEMBER TO CLEAN ROWS AND COLUMNS FOR THE BOUNDARY ELEMENTS !!!
      /*int currBCElNode = 0;
      for(size_t loopA=0;loopA<model->bcElementList.size();loopA++){
        for(size_t loopB=0;loopB<model->bcElementList[loopA]->elementConnections.size();loopB++){
          currBCElNode = model->bcElementList[loopA]->elementConnections[loopB];
          printf("Removed: %d\n",currBCElNode);
          advDiffMat->clearRowAndColumn(currBCElNode);
          advDiffVec->values[currBCElNode] = 0.0;
        }
      }*/

      // LOOP ON BOUNDARY ELEMENTS
      int parentElement = 0;
      for(size_t loopElement=0;loopElement<model->bcElementList.size();loopElement++){

        // Get The Parent Element for this BC Element
        parentElement = model->bcParentElement[loopElement];

        // Assemble Advection-Diffusion Matrix
        model->bcElementList[loopElement]->formWeakBC(model->nodeList,
                                                      model->elementList[parentElement],
                                                      rule,
                                                      (femDoubleVec)model->elDiffusivity[parentElement],
                                                      (femDoubleVec)model->elVelocity[parentElement],
                                                      (femDoubleVec)model->bcElementNormal[loopElement],
                                                      model->bcElementValue[loopElement],
                                                      elMat,elRhs);
        //printf("RHS Local\n");
        //for(int loopC=0;loopC<elRhs.size();loopC++){
        //  printf("RHS %d, %f\n",loopC,elRhs[loopC]);
        //}
        // Assemble both LHS and RHS Contributions
        // Assemble Sparse Matrix
        advDiffMat->assemble(elMat,model->elementList[parentElement]->elementConnections);
        // Assemble Source Vector
        advDiffVec->assemble(elRhs,model->elementList[parentElement]->elementConnections);

      }
    }else{
      // USE STRONG BOUNDARY CONDITIONS
      printf("Assembling Strong Boundary Conditions...\n");
      // Sparse Matrix
      advDiffMat->applyDirichelet(model->diricheletBCNode);
      // RHS Vector
      advDiffVec->applyDirichelet(model->diricheletBCNode,model->diricheletBCValues);
    }

    // PRINT LHS MATRIX
    advDiffMat->writeToFile(string("lhsMatrix.txt"));

    // PRINT RHS VECTOR
    advDiffVec->writeToFile(string("rhsVector.txt"));

    // SOLVE LINEAR SYSTEM OF EQUATIONS
    femVector* solution;
#ifdef USE_ARMADILLO
    solution = solveLinearSystem((femDenseMatrix*)advDiffMat,(femDenseVector*)advDiffVec);
#endif
#ifdef USE_CSPARSE
    solution = solveLinearSystem((femSparseMatrix*)advDiffMat,(femDenseVector*)advDiffVec);
#endif
#ifdef USE_TRILINOS
    solution = solveLinearSystem(model->totNodesInProc,(femTrilinosMatrix*)advDiffMat,(femTrilinosVector*)advDiffVec,nodeDOFs);
#endif

    // ADD SOLUTION TO MODEL RESULTS
    femDoubleVec temp;
    femResult* res = new femResult();
    res->label = string("AdvDiffResult");
    res->type = frNode;
    res->numComponents = 1;
    // Assign to values
    for(int loopA=0;loopA<solution->getSize();loopA++){
      printf("Solution %d %f\n",loopA,solution->getComponent(loopA));
      temp.clear();
      temp.push_back(solution->getComponent(loopA));
      res->values.push_back(temp);
    }
    model->resultList.push_back(res);
}

// ======================
// EVAL DISTANCE FUNCTION
// ======================
void evalDistanceFromPoissonSolution(femVector* femSolution, femModel* model,femDoubleVec& cellDistance){
  femDoubleVec shapeFunction;
  femDoubleMat globShDeriv;
  femDoubleVec elemSol;

  double detJ = 0.0;
  double cellSolution = 0;
  double cellSolDeriv[3];
  cellDistance.clear();
  // Loop through the elements to get the global shape function derivative
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Collect Solution at element nodes
    elemSol.clear();
    for(int loopA=0;loopA<model->elementList[loopElement]->elementConnections.size();loopA++){
      elemSol.push_back(femSolution->getComponent(model->elementList[loopElement]->elementConnections[loopA]));
    }

    // Eval Shape functions and their global derivatives at the cell centers
    model->elementList[loopElement]->evalShapeFunction(model->nodeList,0.0,0.0,0.0,shapeFunction);
    model->elementList[loopElement]->evalGlobalShapeFunctionDerivative(model->nodeList,0.0,0.0,0.0,detJ,globShDeriv);

    // Assemble Cell Solution and its derivatives
    cellSolution = 0;
    cellSolDeriv[0] = 0.0;
    cellSolDeriv[1] = 0.0;
    cellSolDeriv[2] = 0.0;
    for(int loopA=0;loopA<model->elementList[loopElement]->elementConnections.size();loopA++){
      cellSolution    += shapeFunction[loopA] * elemSol[loopA];
      cellSolDeriv[0] += globShDeriv[loopA][0] * elemSol[loopA];
      cellSolDeriv[1] += globShDeriv[loopA][1] * elemSol[loopA];
      cellSolDeriv[2] += globShDeriv[loopA][2] * elemSol[loopA];
    }

    // Assemble into cell distance
    cellDistance.push_back(- sqrt(cellSolDeriv[0] * cellSolDeriv[0] + cellSolDeriv[1] * cellSolDeriv[1] + cellSolDeriv[2] * cellSolDeriv[2])
                           + sqrt(cellSolDeriv[0] * cellSolDeriv[0] + cellSolDeriv[1] * cellSolDeriv[1] + cellSolDeriv[2] * cellSolDeriv[2] + 2.0 * cellSolution));
  }
}

// ======================
// SOLVE POISSON EQUATION
// ======================
void femPoissonSolver::solve(femOption* options, femModel* model){

  // Set Average
  const bool removeAVG = false;
  femDoubleVec cellDistance;

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // One degree of freedom nodes in Poisson Solver
  int nodeDOFs = 1;

  // Init Sparse Matrix and Dense Vector
  femMatrix* poissonMat;
  femVector* poissonVec;
#ifdef USE_ARMADILLO
  poissonMat = new femDenseMatrix(model);
  poissonVec = new femDenseVector((int)model->nodeList.size());
#endif
#ifdef USE_CSPARSE  
  poissonMat = new femSparseMatrix(model);
  poissonVec = new femDenseVector((int)model->nodeList.size());
#endif
#ifdef USE_TRILINOS
  poissonVec = new femTrilinosVector(model->totNodesInProc,model->localToGlobalNodes,nodeDOFs);
  poissonMat = new femTrilinosMatrix(model,nodeDOFs);
#endif  

  // Local Element Matrix
  femDoubleMat elMat;
  elMat.resize(kMaxConnections);
  for(int loopA=0;loopA<kMaxConnections;loopA++){
    elMat[loopA].resize(kMaxConnections);
  }

  // ASSEMBLE MATRIX FROM ALL ELEMENTS
  int precentProgress,percentCounted;
  printf("Assembling Matrix...");
  fflush(stdout);
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    precentProgress = (int)(((double)loopElement/(double)model->elementList.size())*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      printf("%d.",precentProgress);
      fflush(stdout);
    }
    // Gauss Points Loop
    // Assemble Poisson Matrix
    model->elementList[loopElement]->formPoissonMatrix(model->nodeList,rule,(femDoubleVec)model->elDiffusivity[loopElement],elMat);
    // Assemble Sparse Matrix
    poissonMat->assemble(elMat,model->elementList[loopElement]->elementConnections);    
  }
#ifdef USE_TRILINOS
  poissonMat->completeFill();
#endif
  printf("100.OK\n");
  fflush(stdout);

  // Print Time
  //femTime::printToScreen();

  // Local Source Vector
  femDoubleVec elSourceVec;
  elSourceVec.resize(kMaxConnections);
  // ASSEMBLE SOURCE TERM
  printf("Assembling Source...");
  fflush(stdout);
  int currEl = 0;
  double currValue = 0.0;  
  for(size_t loopSource=0;loopSource<model->sourceElement.size();loopSource++){
    precentProgress = (int)(((double)loopSource/(double)model->sourceElement.size())*100);
    if (((precentProgress % 10) == 0)&&((precentProgress / 10) != percentCounted)){
      percentCounted = (precentProgress / 10);
      printf("%d.",precentProgress);
      fflush(stdout);
    }
    // Get Current Element
    currEl = model->sourceElement[loopSource];
    // Eval Source Vector
    model->elementList[currEl]->formPoissonSource(model->nodeList,rule,model->sourceValues[loopSource],elSourceVec);
    // Assemble Source Vector
    poissonVec->assemble(elSourceVec,model->elementList[currEl]->elementConnections);
  }
#ifdef USE_TRILINOS
  poissonVec->GlobalAssemble();
#endif
  printf("100.OK\n");
  fflush(stdout);

  double sourceSum = 0.0;
  for(int loopA=0;loopA<poissonVec->getSize();loopA++){
    sourceSum += poissonVec->getComponent(loopA);
  }
  printf("Source Summation: %f\n",sourceSum);

  // Computing Volume
  double totalVolume = 0.0;
  double currentVolume = 0.0;
  double maxVolume = -std::numeric_limits<double>::max();
  double minVolume = std::numeric_limits<double>::max();
  for(int loopA=0;loopA<model->elementList.size();loopA++){
    currentVolume = model->elementList[loopA]->EvalVolume(0.0,model->nodeList);
    totalVolume += currentVolume;
    if(currentVolume > maxVolume){
      maxVolume = currentVolume;
    }
    if(currentVolume < minVolume){
      minVolume = currentVolume;
    }
  }
  printf("Min Volume: %f, Max Volume: %f\n",minVolume,maxVolume);
  printf("Total Volume: %f\n",totalVolume);


  // ASSEMBLE NEUMANN BOUNDARY CONDITIONS
  printf("Assembling Neumann Vector...\n");
  fflush(stdout);
  femDoubleVec elBCVec;
  int currNode = 0;
  elBCVec.resize(model->nodeList.size());
  for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
    elBCVec[loopA] = 0.0;
  }
  for(size_t loopBC=0;loopBC<model->neumannBCElement.size();loopBC++){
    // Get Current Element
    currEl = model->neumannBCElement[loopBC];
    // Get Value
    currValue = model->neumannBCValues[loopBC];
    // Get Global Face Nodes
    for(size_t loopA=0;loopA<model->neumannBCFaceNodes[loopBC].size();loopA++){
      currNode = model->neumannBCFaceNodes[loopBC][loopA];
      elBCVec[currNode] += (currValue)/(double)model->neumannBCFaceNodes[loopBC].size();
    }
  }

  // Compute Sum of Neumann BC
  double neuBCSum = 0.0;
  for(int loopA=0;loopA<elBCVec.size();loopA++){
    neuBCSum += elBCVec[loopA];
  }
  printf("Neumann BC Summation: %f\n",neuBCSum);

  // Assemble in RHS
  printf("Assembling Neumann RHS...\n");
  fflush(stdout);
  for(int loopA=0;loopA<poissonVec->getSize();loopA++){
    poissonVec->addComponent(loopA,elBCVec[loopA]);
  }
#ifdef USE_TRILINOS
  poissonVec->GlobalAssemble();
#endif

  // APPLY DIRICHELET BOUNDARY CONDITIONS
  if(model->diricheletBCNode.size() > 0){
    printf("Assembling Dirichelet...\n");
    // Sparse Matrix
    poissonMat->applyBlockDirichelet(model->diricheletBCNode,0);
    // RHS Vector
    poissonVec->applyBlockDirichelet(model->diricheletBCNode,model->diricheletBCValues,0);
  }

  // PRINT MATRIX AND VECTOR TO FILE
  //poissonMat->writeToFile(string("pMatrix.txt"));
  //poissonVec->writeToFile(string("pVector.txt"));

  // Check the sum of terms in the matrix
  double sum = 0.0;
  for(int loopA=0;loopA<poissonVec->getSize();loopA++){
    sum += poissonVec->getComponent(loopA);
  }
  printf("RHS Summation: %e\n",sum);

  // SOLVE LINEAR SYSTEM OF EQUATIONS
  femVector* solution;
  printf("Solving Linear System...\n");
  fflush(stdout);
#ifdef USE_ARMADILLO
  solution = solveLinearSystem((femDenseMatrix*)poissonMat,(femDenseVector*)poissonVec);
#endif
#ifdef USE_CSPARSE
  solution = solveLinearSystem((femSparseMatrix*)poissonMat,(femDenseVector*)poissonVec);
#endif
#ifdef USE_TRILINOS
  solution = solveLinearSystem(model->totNodesInProc,(femTrilinosMatrix*)poissonMat,(femTrilinosVector*)poissonVec,1);
#endif

  // EVAL SOLUTION AVERAGE
  double solAVG = 0.0;
  if(model->problemType == ptPPE){
    for(size_t loopA=0;loopA<solution->getSize();loopA++){
      // printf("Solution Comp: %f\n",solution->getComponent(loopA));
      solAVG += solution->getComponent(loopA);
    }
    solAVG = solAVG/(double)solution->getSize();
   printf("Solution Size: %d, Avg: %f\n",solution->getSize(),solAVG);
  }

  // EVAL DISTANCE FUNCTION
  femDoubleVec temp;
  if(model->problemType == ptPoissonDistance){
    // EVALUATE WALL DISTANCE FROM POISSON EQUATION
    evalDistanceFromPoissonSolution(solution,model,cellDistance);
    // ADD CELL DISTANCE SOLUTION
    femResult* cellRes = new femResult();
    cellRes->label = string("cellDistance");
    cellRes->type = frElement;
    cellRes->numComponents = 1;
    // Assign to values
    for(size_t loopA=0;loopA<cellDistance.size();loopA++){
      temp.clear();
      temp.push_back(cellDistance[loopA]);
      cellRes->values.push_back(temp);
    }
    model->resultList.push_back(cellRes);
    // Export Vector To File
    femUtils::writeVectorToFile(string("cellDistance.txt"), cellDistance);
  }

  // ADD SOLUTION TO MODEL RESULTS AND SUBTRACT AVERAGE
  femResult* res = new femResult();
  res->label = string("PoissonResult");
  res->type = frNode;
  res->numComponents = 1;
  // Assign to values
  for(size_t loopA=0;loopA<solution->getSize();loopA++){
    temp.clear();
    temp.push_back(solution->getComponent(loopA) - solAVG);
    res->values.push_back(temp);
  }
  model->resultList.push_back(res);

  // Free Memory
  delete poissonMat;
  delete poissonVec;
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
  for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
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
  femDoubleMat intCoords = rule->getCoords(numNodes,model->elementList[currElement]->dims);

  // Eval Current Shape Derivatives Matrix at the first Gauss Point
  femDoubleMat shapeDeriv;
  double detJ;
  model->elementList[currElement]->evalGlobalShapeFunctionDerivative(model->nodeList,intCoords[0][0],intCoords[0][1],intCoords[0][2],detJ,shapeDeriv);

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

#ifdef USE_TRILINOS
void femSteadyStateAdvectionDiffusionSolver::assembleLHS(femOption* options, femModel* model,Epetra_FECrsMatrix &lhs){
  printf("Assembling Advection-Diffusion LHS...\n");

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Init Sparse Matrix and Dense Vector
  femMatrix* advDiffMat;
  femVector* advDiffVec;
  advDiffMat = new femSparseMatrix(model);
  advDiffVec = new femDenseVector((int)model->nodeList.size());

  // Local Element Matrix
  femDoubleMat elMat;
  femDoubleVec elRhs;
  int lcErr;

  // ASSEMBLE LHS MATRIX
  int totNodes = 0;
  printf("Assembling LHS Matrix...\n");
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Get Total Number of Nodes
    totNodes = model->elementList[loopElement]->elementConnections.size();
    // Assemble Advection-Diffusion Matrix
    model->elementList[loopElement]->formAdvDiffLHS(model->nodeList,
                                                    rule,
                                                    (femDoubleVec)model->elDiffusivity[loopElement],
                                                    (femDoubleVec)model->elVelocity[loopElement],
                                                    ((femAdvectionDiffusionOptions*)options)->advDiffScheme,
                                                    elMat);
    // Fill Trilinos LHS
    // Store the Dense Matrix
    Epetra_SerialDenseMatrix k(totNodes,totNodes);
    for(int loopA=0;loopA<totNodes;loopA++){
      for(int loopB=0;loopB<totNodes;loopB++){
        k(loopA,loopB) = elMat[loopA][loopB];
      }
    }
    Epetra_IntSerialDenseVector indices(totNodes);
    for(int loopA=0;loopA<totNodes;loopA++){
      indices(loopA) = model->elementList[loopElement]->elementConnections[loopA];
    }
    lcErr = lhs.SumIntoGlobalValues(indices,k);
    cout << lcErr << endl;
  }
}
void femSteadyStateAdvectionDiffusionSolver::assembleRHS(femOption* options, femModel* model,Epetra_FEVector &rhs){
  // Local Element Matrix
  femDoubleVec elRhs;

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // ASSEMBLE RHS TERM
  printf("Assembling RHS Vector...\n");
  int currEl = 0;
  int totNodes = 0;
  double currDiff = 0.0;
  double currValue = 0.0;
  for(size_t loopSource=0;loopSource<model->sourceElement.size();loopSource++){
    // Get Current Element
    currEl = model->sourceElement[loopSource];

    // Get Total Nodes
    totNodes = model->elementList[currEl]->elementConnections.size();

    model->sourceValues[loopSource],
    (femDoubleVec)model->elDiffusivity[loopSource],
    (femDoubleVec)model->elVelocity[loopSource],
    ((femAdvectionDiffusionOptions*)options)->advDiffScheme,

    // Eval Source Vector
    model->elementList[currEl]->formAdvDiffRHS(model->nodeList,
                                               rule,
                                               model->sourceValues[loopSource],(femDoubleVec)model->elDiffusivity[loopSource],(femDoubleVec)model->elVelocity[loopSource],
                                               ((femAdvectionDiffusionOptions*)options)->advDiffScheme,
                                               elRhs);
    // Fill Trilinos LHS
    // Store the Dense Matrix
    Epetra_SerialDenseVector k(totNodes);
    for(int loopA=0;loopA<totNodes;loopA++){
      k[loopA] = elRhs[loopA];
    }
    Epetra_IntSerialDenseVector indices(totNodes);
    for(int loopA=0;loopA<totNodes;loopA++){
      indices[loopA] = model->elementList[currEl]->elementConnections[loopA];
    }
    rhs.SumIntoGlobalValues(indices,k);
  }
}
#endif

