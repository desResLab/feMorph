# include "femIncompressibleSolver.h"
# include "femModel.h"
# include "femMatrix.h"
# include "femVector.h"
# include "femSolver.h"
# include "femUtils.h"

#ifdef USE_TRILINOS
# include "trilinos/femTrilinosMatrix.h"
# include "trilinos/femTrilinosVector.h"
#endif

// CONSTRUCTOR
femIncompressibleSolver::femIncompressibleSolver(){
}

// ===================
// PRESCIBE VELOCITIES
// ===================
void prescribeNodeVels(femModel* model, int prescribedVelType,double currTime,femDoubleMat& solution){
  // Rotating Velocities for HW3
  double nodeX = 0.0;
  double nodeY = 0.0;
  if(prescribedVelType == 0){
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
      nodeX = model->nodeList[loopA]->coords[0];
      nodeY = model->nodeList[loopA]->coords[1];
      solution[loopA][0] = cos(M_PI * currTime / 8.0) * sin(2.0 * M_PI * nodeY) * sin(M_PI * nodeX) * sin(M_PI * nodeX);
      solution[loopA][1] = -cos(M_PI * currTime / 8.0) * sin(2.0 * M_PI * nodeX) * sin(M_PI * nodeY) * sin(M_PI * nodeY);
      solution[loopA][2] = 0.0;
    }
  }else if(prescribedVelType == 1){
    // Constant Velocity equal to the prescribed velocity
    // Assume all the Element Velocities are Specified
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
        nodeX = model->nodeList[loopA]->coords[0];
        nodeY = model->nodeList[loopA]->coords[1];
        solution[loopA][0] =  sin(M_PI * nodeX) * cos(M_PI * nodeY);
        solution[loopA][1] = -cos(M_PI * nodeX) * sin(M_PI * nodeY);
        solution[loopA][2] = 0.0;
    }
  }else if(prescribedVelType == 2){
    // Constant Velocity equal to the prescribed velocity
    // Assume all the Element Velocities are Specified
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
      solution[loopA][0] = model->elVelocity[loopA][0];
      solution[loopA][1] = model->elVelocity[loopA][1];
      solution[loopA][2] = model->elVelocity[loopA][2];
    }
  }
}

// ===========================================
// PERFORM ONE INCOMPRESSIBLE NS STEP IN TIME
// ===========================================
void advanceNavierStokes(long loopStep, double currTime,
                         femModel* model,
                         femDoubleMat solution,
                         femDoubleMat solution_Dot,
                         int nodeDOFs,
                         femVector* sol){
  // Set Scheme Flag
  int schemeType = 0;

  // Print Step Info
  printf("Incompressible NS Step %d, Current Time: %f\n",(int)(loopStep+1),currTime);

  // Get Quantities From Options
  double timeStep = model->timeStep;
  double alphaM = model->alphaM;
  double alphaF = model->alphaF;
  double gamma = model->gamma;

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Init Sparse Matrix and Dense Vector
  printf("Assembling Matrix...");
  fflush(stdout);
  femMatrix* nsLHS;
  femVector* nsRHS;
#ifdef USE_ARMADILLO
  nsLHS = new femDenseMatrix(model);
  nsRHS = new femDenseVector((int)model->nodeList.size());
#endif  
#ifdef USE_CSPARSE
  nsLHS = new femSparseMatrix(model);
  nsRHS = new femDenseVector((int)model->nodeList.size());
#endif    
#ifdef USE_TRILINOS  
  nsLHS = new femTrilinosMatrix(model,nodeDOFs);
  nsRHS = new femTrilinosVector(model->totNodesInProc,model->localToGlobalNodes,nodeDOFs);
#endif  
  printf("Done.\n");
  fflush(stdout);

  // INIT LOCAL ELEMENT MATRIX
  femDoubleBlockMat elMat;
  femDoubleBlockVec elRhs;

  // ASSEMBLE LHS MATRIX
  printf("Assembling NS LHS...");
  fflush(stdout);
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Assemble LHS
    model->elementList[loopElement]->formNS_LHS(model->nodeList,
                                                rule,
                                                model->elDensity[loopElement],
                                                model->elViscosity[loopElement],
                                                solution,
                                                schemeType,
                                                nodeDOFs,
                                                timeStep,
                                                alphaM,alphaF,gamma,
                                                elMat);
    // Assemble Sparse Matrix
    nsLHS->blockAssemble(elMat,model->elementList[loopElement]->elementConnections);
  }
  printf("Done.\n");
  fflush(stdout);

  // ASSEMBLE RHS TERM
  // Carefull: No Source Considered!
  printf("Assembling NS RHS...");
  fflush(stdout);
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Eval Source Vector
    model->elementList[loopElement]->formNS_RHS(model->nodeList,
                                                rule,
                                                0.0,
                                                (femDoubleVec)model->elDiffusivity[loopElement],
                                                solution,
                                                solution_Dot,
                                                timeStep,
                                                alphaM,alphaF,gamma,
                                                elRhs);
    // Assemble Source Vector
    nsRHS->blockAssemble(elRhs,model->elementList[loopElement]->elementConnections);
  }
  printf("Done.\n");
  fflush(stdout);

  // Sparse Matrix
  nsLHS->applyDirichelet(model->diricheletBCNode);
  // RHS Vector
  nsRHS->applyDirichelet(model->diricheletBCNode,model->diricheletBCValues);

  // SOLVE LINEAR SYSTEM OF EQUATIONS AND COMPUTE VARIABLE INCREMENT
  printf("Solution...");
  fflush(stdout);
  femSolver linSolver;
  //sol = linSolver.solveLinearSystem((femTrilinosMatrix*)nsLHS,(femTrilinosVector*)nsRHS);
  printf("Done.\n");
  fflush(stdout);
}


// ===========================================
// PERFORM ONE ADVECTON-DIFFUSION STEP IN TIME
// ===========================================
void advanceAdvectionDiffusion(long loopStep, double currTime,
                               femModel* model,
                               femDoubleMat solution,
                               femDoubleMat solution_Dot,
                               int variableID,
                               femVector* sol){

  // Print Step Info
  printf("Advection-Diffusion Step %d, Current Time: %f\n",(int)(loopStep+1),currTime);

  // One degree of freedom in Advection Diffusion
  int nodeDOFs = 1;

  // Get Quantities From Options
  double timeStep = model->timeStep;
  double alphaM = model->alphaM;
  double alphaF = model->alphaF;
  double gamma = model->gamma;

  // Get Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Init Sparse Matrix and Dense Vector
  printf("Assembling Matrix...");
  fflush(stdout);
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
  // One degree of freedom per node in advection diffusion solver
  advDiffVec = new femTrilinosVector(model->totNodesInProc,model->localToGlobalNodes,nodeDOFs);
#endif

  printf("Done.\n");
  fflush(stdout);

  // INIT LOCAL ELEMENT MATRIX
  femDoubleMat elMat;
  femDoubleVec elRhs;

  // ASSEMBLE LHS MATRIX
  printf("Assembling LHS...");
  fflush(stdout);
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Assemble Advection-Diffusion Matrix
    model->elementList[loopElement]->formTransientAdvDiffLHS(model->nodeList,
                                                             rule,
                                                             (femDoubleVec)model->elDiffusivity[loopElement],
                                                             solution,
                                                             timeStep,
                                                             alphaM,alphaF,gamma,
                                                             elMat);
    // Assemble Sparse Matrix
    advDiffMat->assemble(elMat,model->elementList[loopElement]->elementConnections);
  }
  printf("Done.\n");
  fflush(stdout);

  // ASSEMBLE RHS TERM
  // Carefull: No Source Considered!
  printf("Assembling RHS...");
  fflush(stdout);
  for(size_t loopElement=0;loopElement<model->elementList.size();loopElement++){
    // Eval Source Vector
    model->elementList[loopElement]->formTransientAdvDiffRHS(model->nodeList,
                                                             rule,
                                                             0.0,
                                                             (femDoubleVec)model->elDiffusivity[loopElement],
                                                             variableID,
                                                             solution,
                                                             solution_Dot,
                                                             timeStep,
                                                             alphaM,alphaF,gamma,
                                                             elRhs);
    // Assemble Source Vector
    advDiffVec->assemble(elRhs,model->elementList[loopElement]->elementConnections);
  }
  printf("Done.\n");
  fflush(stdout);

    // APPLY DIRICHELET/ESSENTIAL BOUNDARY CONDITIONS

    if(false){
      // USE WEAK BOUNDARY CONDITIONS
      printf("Assembling Weak Boundary Conditions...\n");

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
        // Assemble both LHS and RHS Contributions
        // Assemble Sparse Matrix
        advDiffMat->assemble(elMat,model->elementList[parentElement]->elementConnections);
        // Assemble Source Vector
        advDiffVec->assemble(elRhs,model->elementList[parentElement]->elementConnections);

      }
    }else{
      // Sparse Matrix
      advDiffMat->applyDirichelet(model->diricheletBCNode);
      // RHS Vector
      advDiffVec->applyDirichelet(model->diricheletBCNode,model->diricheletBCValues);
    }

    // PRINT LHS MATRIX
    //advDiffMat->writeToFile(string("lhsMatrix.txt"));

    // PRINT RHS VECTOR
    //advDiffVec->writeToFile(string("rhsVector.txt"));

    // SOLVE LINEAR SYSTEM OF EQUATIONS AND COMPUTE VARIABLE INCREMENT
    printf("Solution...");
    fflush(stdout);
    femSolver linSolver;
#ifdef USE_ARMADILLO
    sol = linSolver.solveLinearSystem((femDenseMatrix*)advDiffMat,(femDenseVector*)advDiffVec);
#endif
#ifdef USE_CSPARSE
    sol = linSolver.solveLinearSystem((femSparseMatrix*)advDiffMat,(femDenseVector*)advDiffVec);
#endif
#ifdef USE_TRILINOS
    sol = linSolver.solveLinearSystem(model->totNodesInProc,(femTrilinosMatrix*)advDiffMat,(femTrilinosVector*)advDiffVec,1);
#endif
    printf("Done.\n");
    fflush(stdout);
}

// SET INITIAL CONDITIONS FOR ALL VARIABLES
void setInitialConditions(femModel* model, femDoubleMat& solution,femDoubleMat& solution_Dot){
  // Allocate and Initialize iniConditions
  int totNodes = model->nodeList.size();
  solution.resize(totNodes);
  solution_Dot.resize(totNodes);
  for(int loopA=0;loopA<totNodes;loopA++){
    solution[loopA].resize(model->maxNodeDofs);
    solution_Dot[loopA].resize(model->maxNodeDofs);
    for(int loopB=0;loopB<model->maxNodeDofs;loopB++){
      solution[loopA][loopB] = 0.0;
      solution_Dot[loopA][loopB] = 0.0;
    }
  }

  // Get Model Initial Conditions
  int currNode = 0;
  int currDof = 0;
  double currVal = 0.0;
  for(size_t loopA=0;loopA<model->iniNodeNumbers.size();loopA++){
    currNode = model->iniNodeNumbers[loopA];
    currDof = model->iniDofNumber[loopA];
    currVal = model->iniDofValue[loopA];
    solution[currNode][currDof] = currVal;
  }

  // Prescribe Initial Velocities
  prescribeNodeVels(model,model->prescribedVelType,0.0,solution);

}

// SOLVE INCOMPRESSIBLE NS WITH ADVECTION STEP
void femIncompressibleSolver::solve(femModel* model){
  // SET UP INITIAL CONDITIONS
  double currentTime = 0.0;

  // Compute total degrees of freedom
  long totNodes = model->nodeList.size();

  // Allocate Solution Vector for this cpu
  // Solution and time derivative and step n
  femDoubleMat solution_n;
  femDoubleMat solution_Dot_n;
  // Solution and time derivative and step n+1
  femDoubleMat solution_n1;
  femDoubleMat solution_Dot_n1;
  // Single timeStep solution for NS
  femVector* timeStepNSSolution;
  // Single timeStep solution for Advection Diffusion
  femVector* timeStepSolution;
  // Quantity for advection diffusion
  int currAdvectedQty = 4;
  // Nodal DOFS for NS Solver
  int nodeDOFs = 4;

  // Temporary Container
  femDoubleVec temp;

  solution_n.resize(totNodes);
  solution_Dot_n.resize(totNodes);
  solution_n1.resize(totNodes);
  solution_Dot_n1.resize(totNodes);
  for(int loopA=0;loopA<totNodes;loopA++){
    // Timestep n
    solution_n[loopA].resize(model->maxNodeDofs);
    solution_Dot_n[loopA].resize(model->maxNodeDofs);
    // Timestep n+1
    solution_n1[loopA].resize(model->maxNodeDofs);
    solution_Dot_n1[loopA].resize(model->maxNodeDofs);
  }

  // SET INITIAL CONDITIONS FOR ALL VARIABLES
  setInitialConditions(model,solution_n,solution_Dot_n);

  // CREATE MODEL RESULTS
  // Create Result for Solution
  femResult* res = new femResult();
  res->label = string("INSSolution");
  res->type = frNode;
  res->numComponents = 1;
  // Assign to values
  for(size_t loopA=0;loopA<solution_n.size();loopA++){
    temp.clear();
    temp.push_back(solution_n[loopA][4]);
    res->values.push_back(temp);
  }
  model->resultList.push_back(res);

  // Create Result for Solution Time Derivative
  femResult* res1 = new femResult();
  res1->label = string("INSSolution_Dot");
  res1->type = frNode;
  res1->numComponents = 1;
  // Assign to values
  for(size_t loopA=0;loopA<solution_Dot_n.size();loopA++){
    temp.clear();
    temp.push_back(solution_Dot_n[loopA][4]);
    res1->values.push_back(temp);
  }
  model->resultList.push_back(res1);

  // Create Result for Velocity
  femResult* res2 = new femResult();
  res2->label = string("velocity");
  res2->type = frNode;
  res2->numComponents = 3;
  // Assign to values
  for(size_t loopA=0;loopA<solution_n.size();loopA++){
    temp.clear();
    temp.push_back(solution_n[loopA][0]);
    temp.push_back(solution_n[loopA][1]);
    temp.push_back(solution_n[loopA][2]);
    res2->values.push_back(temp);
  }
  model->resultList.push_back(res2);

  // Assign to values
  for(size_t loopA=0;loopA<solution_n.size();loopA++){
    model->resultList[0]->values[loopA][0] = solution_n[loopA][4];
    model->resultList[1]->values[loopA][0] = solution_Dot_n[loopA][4];
    model->resultList[2]->values[loopA][0] = solution_n[loopA][0];
    model->resultList[2]->values[loopA][1] = solution_n[loopA][1];
    model->resultList[2]->values[loopA][2] = solution_n[loopA][2];
  }
  // Export Model to Check
  model->ExportToVTKLegacy(string("out_Step_0.vtk"));

  // TIME LOOP
  int currStageType = 0;
  int saveCounter = 0;
  for(long loopTime = 0;loopTime<model->totalSteps+1;loopTime++){
    // Update Current Time
    if(loopTime > 0){
      currentTime += model->timeStep;
    }
    // Update Save Counter
    saveCounter++;
    // LOOP ON STAGES
    for(size_t loopStage = 0;loopStage<model->solStages.size(); loopStage++){

      // Get Stage Type
      currStageType = model->solStages[loopStage];

      // Perform Different Solution based on Stages
      if(currStageType == 0){
        if(model->usePrescribedVelocity){

          // Use prescribed velocities
          printf("Update Velocities...");
          prescribeNodeVels(model,model->prescribedVelType,currentTime,solution_n);
          printf("Done.\n");

        }else{

          // Solve Navier-Stokes Step
          advanceNavierStokes(loopTime,currentTime,
                              model,
                              solution_n,
                              solution_Dot_n,
                              nodeDOFs,
                              timeStepNSSolution);

        }
      }else{

        // Solve Advection-Diffusion Step
        currAdvectedQty = 4;
        advanceAdvectionDiffusion(loopTime,currentTime,
                                  model,
                                  solution_n,
                                  solution_Dot_n,
                                  currAdvectedQty,
                                  timeStepSolution);
      }
    }

    // UPDATE SOLUTION AND TIME DERIVATIVE
    // Update time derivative
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
      solution_Dot_n1[loopA][currAdvectedQty] = timeStepSolution->getComponent(loopA);
    }
    // Update Solution
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
      solution_n1[loopA][currAdvectedQty] = solution_n[loopA][currAdvectedQty] +
                                            model->timeStep * solution_Dot_n[loopA][currAdvectedQty]  +
                                            model->gamma * model->timeStep * (solution_Dot_n1[loopA][currAdvectedQty] - solution_Dot_n[loopA][currAdvectedQty]);
    }
    // Update Old Solutions and time derivatives
    for(size_t loopA=0;loopA<model->nodeList.size();loopA++){
      solution_n[loopA][currAdvectedQty] = solution_n1[loopA][currAdvectedQty];
      solution_Dot_n[loopA][currAdvectedQty] = solution_Dot_n1[loopA][currAdvectedQty];
    }

    // SAVE RESULTS
    if(saveCounter == model->saveEvery){
      // Restore saveCounters
      saveCounter = 0;
      // Assign to values
      for(size_t loopA=0;loopA<solution_n.size();loopA++){
        model->resultList[0]->values[loopA][0] = solution_n[loopA][4];
        model->resultList[1]->values[loopA][0] = solution_Dot_n[loopA][4];
        model->resultList[2]->values[loopA][0] = solution_n[loopA][0];
        model->resultList[2]->values[loopA][1] = solution_n[loopA][1];
        model->resultList[2]->values[loopA][2] = solution_n[loopA][2];
      }
      // Export Model to Check
      model->ExportToVTKLegacy(string("out_Step_" + femUtils::intToStr(loopTime+1) + ".vtk"));
    }
  }
}

