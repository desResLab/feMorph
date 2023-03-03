#include "femGPUFluidSolver.h"

/* Explicit VMS on GPU kernel programs.
 */
const double alpha = 0.58541020;
const double beta = 0.13819660;
const double w[] = {0.25, 0.25, 0.25, 0.25};
const double lN[4][4] = {{0.58541020, 0.13819660, 0.13819660, 0.13819660},
                         {0.13819660, 0.58541020, 0.13819660, 0.13819660},
                         {0.13819660, 0.13819660, 0.58541020, 0.13819660},
                         {0.13819660, 0.13819660, 0.13819660, 0.58541020}};
const double lDN[3][4] = {{-1.0, 1.0, 0.0, 0.0},
                          {-1.0, 0.0, 1.0, 0.0},
                          {-1.0, 0.0, 0.0, 1.0}};
                          
// CONSTRUCTOR
femGPUFluidSolver::femGPUFluidSolver(){
}

/* Initialize the solver, calculate volume, global DN, for each element;
   Assemble the mass matrix.
 */

void initial_assemble(const uint iElm, const femModel *model,
                      femDoubleVec& volumes, femDoubleVec& DNs, femDoubleMat& lumpLHS)
{
    long nNodes = model->nodeList.size();
    long nElms = model->elementList.size();

    long nodeIds[4];
    double nodeCoords[4][3]; // (4,3)

    double jac[3][3];
    double cof[3][3];
    double invJac[3][3];
    double detJ = 0.0;
    double iDetJ = 0.0;

    double value = 0.0;

    // Remember element's nodeIds and nodeCoords.
    for (uint i = 0; i < 4; ++i)
    {
      nodeIds[i] = model->elementList[iElm]->elementConnections[i];

      for (uint j = 0; j < 3; ++j)
      {
        nodeCoords[i][j] = model->nodeList[nodeIds[i]]->coords[j];
      }
    }

    // Calculate jacobian and inverse of jacobian.
    for (uint i = 0; i < 3; ++i)
    {
      for (uint j = 0; j < 3; ++j)
      {
        jac[i][j] = nodeCoords[(j+1)][i] - nodeCoords[0][i];
      }
    }

    // cof[0] = jac[4]*jac[8] - jac[5]*jac[7];
    // cof[1] = jac[5]*jac[6] - jac[3]*jac[8];
    // cof[2] = jac[3]*jac[7] - jac[4]*jac[6];

    // cof[3] = jac[2]*jac[7] - jac[1]*jac[8];
    // cof[4] = jac[0]*jac[8] - jac[2]*jac[6];
    // cof[5] = jac[1]*jac[6] - jac[0]*jac[7];

    // cof[6] = jac[1]*jac[5] - jac[2]*jac[4];
    // cof[7] = jac[2]*jac[3] - jac[0]*jac[5];
    // cof[8] = jac[0]*jac[4] - jac[1]*jac[3];

    // detJ = jac[0]*cof[0] + jac[1]*cof[1] + jac[2]*cof[2];
    // iDetJ = 1.0 / detJ;

    cof[0][0] = jac[1][1]*jac[2][2] - jac[2][1]*jac[1][2];
    cof[0][1] = jac[2][0]*jac[1][2] - jac[1][0]*jac[2][2];
    cof[0][2] = jac[1][0]*jac[2][1] - jac[2][0]*jac[1][1];
    cof[1][0] = jac[2][1]*jac[0][2] - jac[0][1]*jac[2][2];
    cof[1][1] = jac[0][0]*jac[2][2] - jac[2][0]*jac[0][2];
    cof[1][2] = jac[2][0]*jac[0][1] - jac[0][0]*jac[2][1];
    cof[2][0] = jac[0][1]*jac[1][2] - jac[1][1]*jac[0][2];
    cof[2][1] = jac[1][0]*jac[0][2] - jac[0][0]*jac[1][2];
    cof[2][2] = jac[0][0]*jac[1][1] - jac[1][0]*jac[0][1];

    detJ = jac[0][0]*cof[0][0] + jac[0][1]*cof[0][1] + jac[0][2]*cof[0][2];
    iDetJ = 1.0 / detJ;

    // if (iElm == 0)
    // {
    //     dbgDetJ[0] = detJ;
    // }

    for (uint i = 0; i < 3; ++i)
    {
      for (uint j = 0; j < 3; ++j)
      {
        invJac[i][j] = cof[j][i] * iDetJ;
      }
    }

    // 'Assemble' volume and DN.
    volumes[iElm] = detJ / 6.0;

    for (uint i = 0; i < 3; ++i)
    {
      for (uint j = 0; j < 4; ++j)
      {        
        DNs[(i*4+j)*nElms+iElm] = lDN[0][j]*invJac[0][i] + lDN[1][j]*invJac[1][i] + lDN[2][j]*invJac[2][i];
      }
    }

    // Assemble lumped mass matrix, e.g. +volume/4.0 to each vertex.
    for (uint i = 0; i < 4; ++i) // 4 vertices (nodes)
    {
      for (uint j = 0; j < 4; ++j) // 4 d.o.f.s
      {
        lumpLHS[nodeIds[i]][j] = detJ / 24.0;
      }
    }
}


void assemble_RHS(const long iElm, femModel *model,
                  // shape function derivatives in DNs
                  const femDoubleVec& volumes, const femDoubleVec& DNs,
                  // previous u^{n}/p^{n} duP, u^{n-1}/p^{n-1} is preDuP
                  const femDoubleMat& sol_n, const femDoubleMat& sol_prev,
                  // sdus sub-grid velocity for the entire mesh
                  const femDoubleVec& sdus, const femDoubleVec& params,
                  // return rhs
                  femDoubleMat& RHS)
{

    long nNodes = model->nodeList.size();
    long nElms = model->elementList.size();

    long nodeIds[4];
    double Ve;
    double DN[3][4]; // 3*4
    
    double f[4][3];
    double sdu[4][3];
    double hdu[4][3]; // 4*3
    double hp[4];
    double ha[4][3];
    double gradHdu[3][3]; // 3*3

    double wGp;
    double hah[3];
    double hph;
    double sduh[3];
    double fh[3];

    double trGradHdu = 0.0;
    double ahGradHu[3];
    double ahDN;
    double sduhDN;
    double lRes[4][4];

    // parameters
    double nu = params[2];
    double invEpsilon = params[4];

    // Memory clear first.
    for (uint i = 0; i < 4; ++i)
    {
        for (uint j = 0; j < 4; ++j)
        {
            lRes[i][j] = 0.0;
        }
    }

    // Loop on element nodes
    for (uint i = 0; i < 4; ++i)
    {
        // nodeIds[i] = elmNodeIds[i*nElms+iElm];
        nodeIds[i] = model->elementList[iElm]->elementConnections[i];

        // Loop on the three components
        for (uint j = 0; j < 3; ++j)
        {
            // Get subgrid velocities in a matrix format (row - node; column - X,Y,Z)
            sdu[i][j] = sdus[(i*3+j)*nElms+iElm];
            // Get forcing in matrix format (row - node; column - X,Y,Z)
            f[i][j] = model->nodeList[nodeIds[i]]->force[j];
            
            // Update velocity predictor
            hdu[i][j] = 1.5*sol_n[nodeIds[i]][j] - 0.5*sol_prev[nodeIds[i]][j];
            // Update convective velocity with predictor plus subscale
            ha[i][j] = hdu[i][j] + sdu[i][j];
        }
        // Update pressure predictor
        hp[i] = 1.5*sol_n[nodeIds[i]][3] - 0.5*sol_prev[nodeIds[i]][3];
    }

    // Get element volume
    Ve = volumes[iElm];

    // Get shape function derivatives for current element
    // Coordinate loop
    for (uint i = 0; i < 3; ++i)
    {
        // Node loop
        for (uint j = 0; j < 4; ++j)
        {
            DN[i][j] = DNs[(i*4+j)*nElms+iElm];
        }
    }

    // Gradient of velocity predictor and its trace at nodes
    for (uint i = 0; i < 3; ++i) // partial u_x, u_y, u_z
    {
        for (uint j = 0; j < 3; ++j) // partial x, y, z
        {
            gradHdu[i][j] = DN[j][0]*hdu[0][i] + DN[j][1]*hdu[1][i] \
                          + DN[j][2]*hdu[2][i] + DN[j][3]*hdu[3][i];
        }

        trGradHdu += gradHdu[i][i];
    }

    // Loop over the Gauss integration points
    for (uint iGp = 0; iGp < 4; ++iGp)
    {
        // Get Gauss point weight  
        wGp = w[iGp] * Ve;

        // Calculate values at Gaussian point iGp.
        for (uint i = 0; i < 3; ++i)
        {
            // Advective velocity at Gauss point
            hah[i] = 0.0;
            // Subgrid scale velocity
            sduh[i] = 0.0;
            // Forcing at the Gauss point
            fh[i] = 0.0;

            for (uint j = 0; j < 4; ++j)
            {
                hah[i] += ha[j][i] * lN[iGp][j];
                sduh[i] += sdu[j][i] * lN[iGp][j];
                fh[i] += f[j][i] * lN[iGp][j];
            }
        }
        
        // Pressure predictor at the Gauss point
        hph = 0.0;
        for (uint i = 0; i < 4; ++i)
        {
            hph += hp[i] * lN[iGp][i];
        }

        // Compute a*grad predictor u
        for (uint i = 0; i < 3; ++i)
        {
            ahGradHu[i] = 0.0;

            for (uint j = 0; j < 3; ++j)
            {
                ahGradHu[i] += gradHdu[i][j] * hah[j];
            }
        }


        // Assemble for each point.
        for (uint a = 0; a < 4; ++a)
        {
            // Calculate temporary value first.
            ahDN = 0.0;
            sduhDN = 0.0;
            for (uint i = 0; i < 3; ++i)
            {
                ahDN += hah[i] * DN[i][a];
                sduhDN += sduh[i] * DN[i][a];
            }
            
            // Assemble first 3 d.o.f.s.
            for (uint i = 0; i < 3; ++i)
            {
                for (uint j = 0; j < 3; ++j)
                {
                    // Assemble viscous term
                    lRes[a][i] += wGp*nu*gradHdu[i][j]*DN[j][a];
                }

                lRes[a][i] += wGp*(ahGradHu[i]*lN[iGp][a] \
                            - hph*DN[i][a] - sduh[i]*ahDN \
                            - fh[i]*lN[iGp][a]);
                
                // lRes[a][i] += wGp*ahGradHu[i]*lN[iGp][a];
            }

            // Assemble last d.o.f. for pressure.
            // Should sduhDN be just DN????
            lRes[a][3] += wGp*(trGradHdu*lN[iGp][a] - sduhDN)*invEpsilon;
        }
    }

    // Assemble to global RHS.
    for (uint a = 0; a < 4; ++a)
    {
        for (uint i = 0; i < 4; ++i)
        {
            RHS[nodeIds[a]][i] += lRes[a][i];
        }
    }
}

// APPLY DIRICHLET BC
void apply_drchBC(femModel* model,femDoubleMat& sol)
{
  for (ulint loopA = 0; loopA < model->diricheletBCNode.size(); loopA++)
  {
    for (ulint loopB = 0; loopB < 3; loopB++)
    {
      sol[model->diricheletBCNode[loopA]][loopB] = model->diricheletBCValues[loopA][loopB];
    }
  }
}

// SOLVE INCOMPRESSIBLE NS WITH ADVECTION STEP
void femGPUFluidSolver::solve(femModel* model){

  // SET UP INITIAL CONDITIONS
  double currentTime = 0.0;

  // Compute total degrees of freedom
  long totNodes = model->nodeList.size();
  long totElements = model->elementList.size();

  // Declare main variables - Set to zero
  femDoubleMat sol_n;
  femDoubleMat sol_prev;
  femDoubleMat rhs;
  
  // SET PRESCRIBED VELOCITIES
  model->prescribeNodeVels(currentTime,sol_n);

  // ASSMBLE ELEMENT QTY
  femDoubleVec volumes(totElements,0.0);
  femDoubleVec DNs(totNodes*model->maxNodeDofs,0.0);
  femDoubleVec sdus(totNodes*model->maxNodeDofs,0.0);
  femDoubleVec params(totNodes*model->maxNodeDofs,0.0);
  femDoubleMat lumpLHS;

  // Assemble initial quantities
  for(uint loopElement=0;loopElement<totElements;loopElement++){
    initial_assemble(loopElement, model,
                     volumes,DNs,lumpLHS);

   // NEED TO UPDATE THE SUBSCALE VELOCITY AND PRESSURE!!!!
  }

  // =========
  // TIME LOOP
  // =========
  long saveCounter = 0;
  for(uint loopTime=0;loopTime<model->totalSteps;loopTime++){
      
    // Set RHS to zero
    for(ulint loopA=0;loopA < totNodes;loopA++){
      for(ulint loopB=0;loopB < 4;loopB++){
        rhs[loopA][loopB] = 0.0;
      }
    }

    // Assemble element contribution in global RHS vector
    for(ulint loopElement=0; loopElement < totElements; loopElement++){

      assemble_RHS(loopElement,model,
                  volumes,DNs,
                  sol_n,sol_prev,
                  sdus,params,
                  rhs);

    }

    // Update main variables
    for(ulint loopNode=0;loopNode<totNodes;loopNode++){
      for(ulint loopDof=0;loopDof<3;loopDof++){
        sol_n[loopNode][loopDof] = sol_prev[loopNode][loopDof] - model->timeStep * rhs[loopNode][loopDof] / lumpLHS[loopNode][loopDof];        
      }
    }

    // APPLY DIRICHLET BC
    apply_drchBC(model,sol_n);

    // Update previous solution
    for(ulint loopNode=0;loopNode<totNodes;loopNode++){
      for(ulint loopDof=0;loopDof<3;loopDof++){
        sol_prev[loopNode][loopDof] = sol_n[loopNode][loopDof];
      }
    }

    // SAVE RESULTS
    if(saveCounter == model->saveEvery){
      // Restore saveCounters
      saveCounter = 0;
      // Assign to values
      for(size_t loopA=0;loopA<sol_n.size();loopA++){
        model->resultList[0]->values[loopA][0] = sol_n[loopA][4];
        model->resultList[1]->values[loopA][0] = sol_n[loopA][0];
        model->resultList[1]->values[loopA][1] = sol_n[loopA][1];
        model->resultList[1]->values[loopA][2] = sol_n[loopA][2];
      }
      // Export Model to Check
      model->ExportToVTKLegacy(string("out_Step_" + femUtils::intToStr(loopTime+1) + ".vtk"));
    }
    // Update the solution vector 

    // Update counter
    saveCounter++;
    // Update current time
    currentTime += model->timeStep;
  }
}


//   // CREATE MODEL RESULTS
//   // Create Result for Solution
//   femResult* res = new femResult();
//   res->label = string("INSSolution");
//   res->type = frNode;
//   res->numComponents = 1;
//   // Assign to values
//   for(size_t loopA=0;loopA<solution_n.size();loopA++){
//     temp.clear();
//     temp.push_back(solution_n[loopA][4]);
//     res->values.push_back(temp);
//   }
//   model->resultList.push_back(res);




