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

// Compute element residual projection
femDoubleMat eval_ortho_proj(const femIntVec& conns, const femDoubleMat& rhs, const double tau_v, const double Ve){
  femDoubleMat ortho_rhs;
  femUtils::matZeros(ortho_rhs,4,3);
  double integral_num = 0.0;
  double integral_den = 0.0;
  double wGp = 0.0;
  double rhs_gauss = 0.0;
  double rhs_ln = 0.0;

  // loop through the nodes (shape function)
  for(ulint loopA=0;loopA<4;loopA++){
    // loop over the global coordinates X,Y,Z
    for(ulint loopB=0;loopB<3;loopB++){

      integral_num = 0.0;
      integral_den = 0.0;
      // Loop over the gauss point
      for(ulint iGp=0;iGp<4;iGp++){

        // Get Gauss point weight  
        wGp = w[iGp] * Ve;

        // Get rhs*lN and lN at the current gauss point
        rhs_gauss = 0.0;
        rhs_ln = 0.0;
        for(ulint i=0;i<4;i++){
          rhs_gauss += rhs[conns[i]][loopB]*lN[i][loopA];
          rhs_ln += lN[i][loopA];
        }

        // Sum gauss point contribution to integral
        integral_num += tau_v*rhs_gauss*rhs_ln*wGp;
        integral_den += rhs_ln*wGp;

      }
      // Assign value to orthogonal 
      if(fabs(integral_den) > kMathZero){
        ortho_rhs[loopA][loopB] = integral_num/integral_den;
      }else{
        throw femException("Zero integral_den in eval_ortho_proj.\n");
      }
    }
  }
  return ortho_rhs;
}

//Initialize the solver, calculate volume, global DN, for each element;
//Assemble the mass matrix.
void initial_assemble(const uint iElm, const femModel *model,
                      femDoubleVec& volumes, femDoubleVec& DNs, femDoubleMat& lumpLHS)
{
    long nNodes = model->nodeList.size();
    long nElms = model->elementList.size();

    long nodeIds[4];
    double nodeCoords[4][3];

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
                  femDoubleMat& sdu,
                  // Return convective velocity
                  femDoubleMat& ha,
                  // return rhs
                  femDoubleMat& RHS)
{

    long nNodes = model->nodeList.size();
    long nElms = model->elementList.size();

    long   nodeIds[4];
    double Ve;
    double DN[3][4]; // 3*4
    
    double f[4][3];
    // double sdu[4][3];
    double hdu[4][3]; // 4*3
    double hp[4];
    // double ha[4][3];
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
    double nu = params[1];
    double invEpsilon = params[2];

    // Memory clear first.
    for (uint i = 0; i < 4; ++i)
    {
        for (uint j = 0; j < 4; ++j)
        {
            lRes[i][j] = 0.0;
        }
    }

    // UPDATE NODAL VARIABLES BEFORE GAUSS POINT LOOP
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

// SOLVE INCOMPRESSIBLE NS WITH ADVECTION STEP
void femGPUFluidSolver::solve(femModel* model){

  // SET UP INITIAL CONDITIONS
  double currentTime = 0.0;

  // Compute total degrees of freedom
  long totNodes = model->nodeList.size();
  long totElements = model->elementList.size();

  // Declare main variables - Set to zero
  femDoubleMat sol_n;
  femUtils::matZeros(sol_n,totNodes,model->maxNodeDofs);
  femDoubleMat sol_prev;
  femUtils::matZeros(sol_prev,totNodes,model->maxNodeDofs);
  femDoubleMat rhs;
  femUtils::matZeros(rhs,totNodes,model->maxNodeDofs);
  // ELEMENT-BASED SUB-GRID SCALE VELOCITY INITIALIZED TO ZERO
  femDoubleVec sdus(totElements*4*3,0.0);

  // SET INITIAL INFLOW VELOCITIES
  model->setNodeVelocity(currentTime,sol_n);

  // SET WALL BOUNDARY CONDITIONS
  // Avoids non-zero velocities at the wall
  model->setDirichletBC(sol_n);

  // COPY PREVIOUS SOLUTION
  // Update previous solution
  for(ulint loopNode=0;loopNode<totNodes;loopNode++){
    for(ulint loopDof=0;loopDof<model->maxNodeDofs;loopDof++){
      sol_prev[loopNode][loopDof] = sol_n[loopNode][loopDof];
    }
  }

  // CREATE MODEL RESULTS
  // Create Result for Solution
  femResult* res = new femResult();
  femDoubleVec temp;
  res->label = string("velocity");
  res->type = frNode;
  res->numComponents = 3;
  // Assign to values
  for(size_t loopA=0;loopA<sol_n.size();loopA++){
    temp.clear();
    temp.push_back(sol_n[loopA][0]);
    temp.push_back(sol_n[loopA][1]);
    temp.push_back(sol_n[loopA][2]);
    res->values.push_back(temp);
  }
  model->resultList.push_back(res);
  res = new femResult();
  res->label = string("pressure");
  res->type = frNode;
  res->numComponents = 1;
  // Assign to values
  for(size_t loopA=0;loopA<sol_n.size();loopA++){
    temp.clear();
    temp.push_back(sol_n[loopA][3]);
    res->values.push_back(temp);
  }
  model->resultList.push_back(res);
  model->ExportToVTKLegacy(string("out_Step_" + femUtils::intToStr(0) + ".vtk"),false);

  // ASSMBLE ELEMENT QTY
  femDoubleVec volumes(totElements,0.0);
  femDoubleVec DNs(totElements*model->maxNodeDofs*3,0.0);
  femDoubleMat lumpLHS;
  femUtils::matZeros(lumpLHS,totNodes,model->maxNodeDofs);

  // Assemble initial quantities
  for(uint loopElement=0;loopElement<totElements;loopElement++){
    initial_assemble(loopElement, model,
                     volumes,DNs,lumpLHS);
  }

  // WRITE CFL MESSAGE
  // Compute minimum element length
  double min_el_h = std::numeric_limits<double>::max();
  for(uint loopElement=0;loopElement<totElements;loopElement++){
    if(pow(6.0*volumes[loopElement],1.0/3.0) < min_el_h){
      min_el_h = pow(6.0*volumes[loopElement],1.0/3.0);
    }
  }
  // Get maximum velocity 
  double max_ini_el_u = 0.0;
  for(uint loopNode=0;loopNode<totNodes;loopNode++){
    if(sqrt(sol_n[loopNode][0]*sol_n[loopNode][0] + sol_n[loopNode][1]*sol_n[loopNode][1] + sol_n[loopNode][2]*sol_n[loopNode][2]) > max_ini_el_u){
      max_ini_el_u = sqrt(sol_n[loopNode][0]*sol_n[loopNode][0] + sol_n[loopNode][1]*sol_n[loopNode][1] + sol_n[loopNode][2]*sol_n[loopNode][2]);
    }
  }

  printf("--- CFL Condition\n");
  printf("Minium element length: %e\n",min_el_h);
  printf("Maximum initial velocity : %e\n",max_ini_el_u);
  printf("CFL number based on assigned initial velocity: %.3f\n",max_ini_el_u*model->timeStep/min_el_h);
  printf("\n");

  // VARIABLES FOR TIME LOOP
  long saveCounter = 0;
  femDoubleMat ha;
  femUtils::matZeros(ha,4,3);
  double temp_double = 0.0;
  double max_conv_a = 0.0;
  femDoubleMat sdu;
  femUtils::matZeros(sdu,4,3);
  // Viscosity
  double nu = model->vmsProps[1];
  // characteristic length
  double el_h = 0.0;
  // Stabilization constants
  double c1 = 4.0;
  double c2 = 2.0;
  double tau1 = 0.0;
  double tau2 = 0.0;
  double tau_v = 0.0;
  double tau_p = 0.0;
  femDoubleMat ortho_rhs;
  // Store element velocity and pressure residual at element centroid
  femDoubleMat el_u_res;
  femUtils::matZeros(el_u_res,totElements,3);
  femDoubleVec el_p_res(totElements,0.0);
  // Centroid velocity and pressure residual norms
  double el_u_res_norm = 0.0;
  double el_p_res_norm = 0.0;
  // Store subgid velocity at element centroid
  femDoubleMat el_u_sdus;
  femUtils::matZeros(el_u_sdus,totElements,3);
  double curr_u_res = 0.0;
  long int curr_node = 0;
  femIntVec conns;

  // =========
  // TIME LOOP
  // =========
  for(uint loopTime=0;loopTime<model->totalSteps;loopTime++){

    // Update current time
    currentTime += model->timeStep;
      
    // Set RHS to zero
    for(ulint loopA=0;loopA<totNodes;loopA++){
      for(ulint loopB=0;loopB<model->maxNodeDofs;loopB++){
        rhs[loopA][loopB] = 0.0;
      }
    }
    
    // INITIALIZE ELEMENT RESIDUALS FOR PLOTTING
    for(ulint loopA=0;loopA<totElements;loopA++){
      el_p_res[loopA] = 0.0;
      for(ulint loopB=0;loopB<3;loopB++){
        el_u_res[loopA][loopB] = 0.0;
        el_u_sdus[loopA][loopB] = 0.0;
      }
    }

    // Assemble element contribution in global RHS vector
    for(ulint loopElement=0;loopElement<totElements;loopElement++){

      // ASSEMLE ELEMENT RESIDUAL
      assemble_RHS(loopElement,model,
                   volumes,DNs,
                   sol_n,sol_prev,
                   sdus,model->vmsProps,
                   // returns
                   sdu,ha,rhs);       

      // COMPUTE RESIDUALS AT ELEMENT CENTROID
      for(ulint loopA=0;loopA<4;loopA++){
        curr_node = model->elementList[loopElement]->elementConnections[loopA];
        for(ulint loopB=0;loopB<3;loopB++){
          el_u_res[loopElement][loopB] += rhs[curr_node][loopB]*0.25;
        }
        el_p_res[loopElement] += rhs[curr_node][3]*0.25;
      }

      // COMPUTE SUBGRID VELOCITY AT ELEMENT CENTROID
      for(ulint loopA=0;loopA<4;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          el_u_sdus[loopElement][loopB] += sdu[loopA][loopB]*0.25;
        }
      }

      // COMPUTE STABILIZATION CONSTANTS FOR CURRENT ELEMENT
      // Max velocity module
      max_conv_a = 0.0;
      for(ulint loopA=0;loopA<4;loopA++){
        temp_double = sqrt(ha[loopA][0]*ha[loopA][0] + ha[loopA][1]*ha[loopA][1] + ha[loopA][2]*ha[loopA][2]);
        if(temp_double > max_conv_a){
          max_conv_a = temp_double;
        }
      }
      
      // Compute characteristic length
      el_h = pow(6.0*volumes[loopElement],1.0/3.0); // Cubit root of 6 times the volume
      // Compute tau1
      tau1 = 1.0/(((c1*nu)/(el_h*el_h)) + ((c2*max_conv_a)/(el_h)));
      tau2 = ((el_h*el_h)/(c1*tau1));
      tau_v = 1.0/((1.0/model->timeStep) + (1.0/tau1));
      tau_p = 1.0/((1.0/model->timeStep) + (1.0/tau2));

      // COMPUTE RESIDUAL PROJECTION
      conns.clear();
      for(ulint loopA=0;loopA<4;loopA++){
        conns.push_back(model->elementList[loopElement]->elementConnections[loopA]);
      }
      ortho_rhs = eval_ortho_proj(conns,rhs,tau_v,volumes[loopElement]);

      // REMOVE PROJECTION FROM ORIGINAL RESIDUAL
      for(ulint loopA=0;loopA<4;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          ortho_rhs[loopA][loopB] = rhs[conns[loopA]][loopB] - ortho_rhs[loopA][loopB];
        }
      }

      // UPDATE SUBGRID VELOCITIES
      // Loop on the element nodes
      for(ulint loopA=0;loopA<4;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          sdus[(loopA*3+loopB)*totElements+loopElement] = (tau_v/model->timeStep)*sdus[(loopA*3+loopB)*totElements+loopElement] - ortho_rhs[loopA][loopB];
        }
      } 
    }

    // COMPUTE MAX RESIDUAL NORMS
    el_u_res_norm = 0.0;
    el_p_res_norm = 0.0;
    for(ulint loopElement=0;loopElement<totElements;loopElement++){
      curr_u_res = sqrt(el_u_res[loopElement][0]*el_u_res[loopElement][0] + el_u_res[loopElement][1]*el_u_res[loopElement][1] + el_u_res[loopElement][2]*el_u_res[loopElement][2]);
      if(curr_u_res > el_u_res_norm){
        el_u_res_norm = curr_u_res;
      }
      if(fabs(el_p_res[loopElement]) > el_p_res_norm){
        el_p_res_norm = fabs(el_p_res[loopElement]);
      }
    }

    // UPDATE MAIN VARIABLES
    for(ulint loopNode=0;loopNode<totNodes;loopNode++){
      for(ulint loopDof=0;loopDof<model->maxNodeDofs;loopDof++){
        sol_n[loopNode][loopDof] = sol_prev[loopNode][loopDof] - model->timeStep * rhs[loopNode][loopDof] / lumpLHS[loopNode][loopDof];        
      }
    }

    // APPLY INLET VELOCITIES
    model->setNodeVelocity(currentTime,sol_n);
    // SET DIRICHLET BC
    model->setDirichletBC(sol_n);

    // WRITE MESSAGES
    if(loopTime == 0){
      printf("STARTING TIME LOOP\n");
      printf("%10s %20s %20s %20s\n","TIME LOOP","TIME [s]","MAX RES NORM U","MAX RES NORM P");
    }
    printf("%10d %20.8f %20.3e %20.3e\n",loopTime+1,currentTime,el_u_res_norm,el_p_res_norm);

    // Update previous solution
    for(ulint loopNode=0;loopNode<totNodes;loopNode++){
      for(ulint loopDof=0;loopDof<model->maxNodeDofs;loopDof++){
        sol_prev[loopNode][loopDof] = sol_n[loopNode][loopDof];
      }
    }

    // SAVE RESULTS
    if(saveCounter == model->saveEvery){      
      // Restore saveCounters
      saveCounter = 0;
      // Reset result list
      model->resultList.clear();
      // SAVE VELOCITY
      res = new femResult();
      res->label = string("velocity");
      res->type = frNode;
      res->numComponents = 3;
      // Assign to values
      for(size_t loopA=0;loopA<sol_n.size();loopA++){
        temp.clear();
        temp.push_back(sol_n[loopA][0]);
        temp.push_back(sol_n[loopA][1]);
        temp.push_back(sol_n[loopA][2]);
        res->values.push_back(temp);
      }
      model->resultList.push_back(res);
      // SAVE PRESSURE
      res = new femResult();
      res->label = string("pressure");
      res->type = frNode;
      res->numComponents = 1;
      // Assign to values
      for(size_t loopA=0;loopA<sol_n.size();loopA++){
        temp.clear();
        temp.push_back(sol_n[loopA][3]);
        res->values.push_back(temp);
      }
      model->resultList.push_back(res);
      // SAVE VELOCITY RESIDUAL 
      res = new femResult();
      res->label = string("res_v");
      res->type = frElement;
      res->numComponents = 3;
      // Assign to values
      for(size_t loopA=0;loopA<totElements;loopA++){
        temp.clear();
        temp.push_back(el_u_res[loopA][0]);
        temp.push_back(el_u_res[loopA][1]);
        temp.push_back(el_u_res[loopA][2]);
        res->values.push_back(temp);
      }
      model->resultList.push_back(res);
      // SAVE PRESSURE RESIDUAL 
      res = new femResult();
      res->label = string("res_p");
      res->type = frElement;
      res->numComponents = 1;
      // Assign to values
      for(size_t loopA=0;loopA<totElements;loopA++){
        temp.clear();
        temp.push_back(el_p_res[loopA]);
        res->values.push_back(temp);
      }
      model->resultList.push_back(res);
      // SAVE SUBGRID VELOCITY
      res = new femResult();
      res->label = string("subgrid_u");
      res->type = frElement;
      res->numComponents = 3;
      // Assign to values
      for(size_t loopA=0;loopA<totElements;loopA++){
        temp.clear();
        temp.push_back(el_u_sdus[loopA][0]);
        temp.push_back(el_u_sdus[loopA][1]);
        temp.push_back(el_u_sdus[loopA][2]);
        res->values.push_back(temp);
      }
      model->resultList.push_back(res);
      // EXPORT MODEL
      model->ExportToVTKLegacy(string("out_Step_" + femUtils::intToStr(loopTime) + ".vtk"),false);
    }

    // Update counter
    saveCounter++;
  }
}


