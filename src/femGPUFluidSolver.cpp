#include "femGPUFluidSolver.h"

// Explicit VMS Solver
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

// Initialize the solver, calculate volume, global DN, for each element;
// Assemble the mass matrix.
void initial_assemble(const uint iElm, const femModel *model,
                      femDoubleVec& char_length,femDoubleVec& volumes, femDoubleVec& DNs,
                      femDoubleVec& nd_mass){

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
  for (uint i = 0; i < 4; ++i){
    nodeIds[i] = model->elementList[iElm]->elementConnections[i];
    for (uint j = 0; j < 3; ++j){
      nodeCoords[i][j] = model->nodeList[nodeIds[i]]->coords[j];
    }
  }

  // Calculate jacobian and inverse of jacobian.
  for (uint i = 0; i < 3; ++i){
    for (uint j = 0; j < 3; ++j){
      jac[i][j] = nodeCoords[(j+1)][i] - nodeCoords[0][i];
    }
  }

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

  for (uint i = 0; i < 3; ++i){
    for (uint j = 0; j < 3; ++j){
      invJac[i][j] = cof[j][i] * iDetJ;
    }
  }

  // 'Assemble' volume and DN.
  volumes[iElm] = detJ / 6.0;

  // Compute characteristic length - CAN IMPROVE WITH CIRC-DIAM!!!
  char_length[iElm]= pow(6.0*volumes[iElm],1.0/3.0); // Cubic root of 6 times the volume

  // Shape function derivatives
  for (uint i = 0; i < 3; ++i){
    for (uint j = 0; j < 4; ++j){
      DNs[(i*4+j)*nElms+iElm] = lDN[0][j]*invJac[0][i] \
                              + lDN[1][j]*invJac[1][i] \
                              + lDN[2][j]*invJac[2][i];
    }
  }

  for (uint i = 0; i < 4; ++i){
    nd_mass[nodeIds[i]] += volumes[iElm]/4.0;
  }
  
}

void assemble_RHS(const long iElm, const long loopTime, femModel *model,
                  // shape function derivatives in DNs
                  const femDoubleVec& volumes, const femDoubleVec& DNs,
                  // previous u^{n}/p^{n} duP, u^{n-1}/p^{n-1} is preDuP
                  const femDoubleMat& sol_n, const femDoubleMat& sol_prev,
                  // sdus sub-grid velocity for the entire mesh
                  const femDoubleVec& sdus, 
                  // sdus sub-grid pressure for the entire mesh
                  const femDoubleVec& sgps, 
                  const femDoubleVec& params,
                  // Return local element-based copy of subgrid-scale velocities
                  femDoubleMat& sdu,
                  // Return convective velocity
                  femDoubleMat& ha,
                  // return rhs
                  femDoubleMat& RHS){

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
  double hah[4][3];
  double hph;
  double sduh[3];
  double fh[3];

  double trGradHdu = 0.0;
  double ahGradHu[3];
  double ahDN;
  double sduhDN;
  double lRes[4][4];

  femDoubleMat u_diff_term;
  femUtils::matZeros(u_diff_term,4,3);
  femDoubleMat u_conv_term;
  femUtils::matZeros(u_conv_term,4,3);
  femDoubleMat u_pressgrad_term;
  femUtils::matZeros(u_pressgrad_term,4,3);
  femDoubleMat u_conv_sg_term;
  femUtils::matZeros(u_conv_sg_term,4,3);
  femDoubleMat u_forcing_term;
  femUtils::matZeros(u_forcing_term,4,3);
  femDoubleMat p_sg_mom_term;
  femUtils::matZeros(p_sg_mom_term,4,3);
  femDoubleVec p_grad_term(4,0.0);
  femDoubleVec p_sg_term(4,0.0);
  
  double sgp = 0.0;

  double u_diff_norm      = 0.0;
  double u_conv_norm      = 0.0;
  double u_pressgrad_norm = 0.0;
  double u_conv_sg_norm   = 0.0;
  double u_forcing_norm   = 0.0;
  double p_grad_norm      = 0.0;
  double p_sg_norm        = 0.0;

  // Set to zero 
  for (uint i = 0; i < 4; i++){
    p_grad_term[i] = 0.0;
    p_sg_term[i] = 0.0;

    for (uint j = 0; j < 3; j++){

      u_diff_term[i][j] = 0.0;
      u_conv_term[i][j] = 0.0;
      u_pressgrad_term[i][j] = 0.0;
      u_conv_sg_term[i][j] = 0.0;
      u_forcing_term[i][j] = 0.0;
      p_sg_mom_term[i][j] = 0.0;
    }
  }
  // parameters
  // Kinematic viscosity: dynamic viscosity/density
  double nu = params[1]/params[0];
  double invEpsilon = params[2];
  
  // Memory clear first.
  for (uint i = 0; i < 4; ++i){ // node
    for (uint j = 0; j < 4; ++j){ // d.o.f.        
      lRes[i][j] = 0.0;
    }
  }

  // Remember element's nodeIds.
  for (uint i = 0; i < 4; ++i){
      nodeIds[i] = model->elementList[iElm]->elementConnections[i];

      // Calculate initial values.
      for (uint j = 0; j < 3; ++j){
          
          // Local copy of the forcing term
          f[i][j] = model->nodeList[nodeIds[i]]->force[j];
          
          // Compute velocity predictor at nodes
          hdu[i][j] = 1.5*sol_n[nodeIds[i]][j] - 0.5*sol_prev[nodeIds[i]][j];
          
          // Value of subgrid velocity at the Gauss points
          sdu[i][j] = sdus[(i*3+j)*nElms+iElm];
      }
      // Compute pressure predictor at nodes
      hp[i] = 1.5*sol_n[nodeIds[i]][3] - 0.5*sol_prev[nodeIds[i]][3];
  }

  // Get element volume
  Ve = volumes[iElm];

  // Local copy of the global SF derivatives at the Gauss points
  for (uint i = 0; i < 3; ++i){
    for (uint j = 0; j < 4; ++j){
      DN[i][j] = DNs[(i*4+j)*nElms+iElm];
    }
  }

  // Calculate Velocity Gradient at the Gauss points - Is a CONSTANT for P1P1
  for (uint i = 0; i < 3; ++i){ // partial u_x, u_y, u_z
    for (uint j = 0; j < 3; ++j){ // partial x, y, z      
          gradHdu[i][j] = DN[j][0]*hdu[0][i] \
                        + DN[j][1]*hdu[1][i] \
                        + DN[j][2]*hdu[2][i] \
                        + DN[j][3]*hdu[3][i];
      }
      // Trace at the Gauss Points
      trGradHdu += gradHdu[i][i];
  }

  // Loop through Gaussian points to do numerical integration.
  for (uint iGp = 0; iGp < 4; ++iGp){
      
    // Get Gauss point weight
    wGp = w[iGp] * Ve;

    // Get sub-grid pressure at current Gauss Point
    sgp = sgps[(iGp)*nElms + iElm];

    // Forcing and convective velocity at iGp.
    for (uint i = 0; i < 3; ++i){

      hah[iGp][i] = 0.0;
      fh[i] = 0.0;

      for (uint j = 0; j < 4; ++j){
        // From velocity predictor
        hah[iGp][i] += hdu[j][i] * lN[iGp][j];
        // From nodal forcing
        fh[i] += f[j][i] * lN[iGp][j];
      }
      // Add subgrid velocity - it is already at the Gauss points
      hah[iGp][i] += sdu[iGp][i];
    }

    // Pressure Predictor at the iGp
    hph = 0.0;
    for (uint i = 0; i < 4; ++i){
      hph += hp[i] * lN[iGp][i];
    }

    // Calculate medium values.
    for (uint i = 0; i < 3; ++i){
      ahGradHu[i] = 0.0;
      for (uint j = 0; j < 3; ++j){
        // OK, this is consistent with the componentwise NS
        ahGradHu[i] += gradHdu[i][j] * hah[iGp][j];
      }
    }

    // Node Loop
    for (uint a = 0; a < 4; ++a){
      // Calculate temporary value first.
      ahDN = 0.0;
      sduhDN = 0.0;
      for (uint i = 0; i < 3; ++i){
        ahDN += hah[iGp][i] * DN[i][a];
        sduhDN += sdu[iGp][i] * DN[i][a];
      }
          
      // Assemble first 3 d.o.f.s.
      for (uint i = 0; i < 3; ++i){
        // Diffusion Term
        for (uint j = 0; j < 3; ++j){
         u_diff_term[a][i] += wGp*nu*(gradHdu[i][j])*DN[j][a];
         lRes[a][i] += wGp*nu*(gradHdu[i][j])*DN[j][a];
        }
        // Convection Term
        u_conv_term[a][i] += wGp*ahGradHu[i]*lN[iGp][a];
        lRes[a][i] +=   wGp*ahGradHu[i]*lN[iGp][a];
      
        // Pressure Gradient Term
        u_pressgrad_term[a][i] += - wGp*hph*DN[i][a];
        lRes[a][i] += - wGp*hph*DN[i][a];

        // SG Convection Term
        u_conv_sg_term[a][i] += - wGp*sdu[iGp][i]*ahDN;
        lRes[a][i] += - wGp*sdu[iGp][i]*ahDN;

        // SG pressure term
        p_sg_mom_term[a][i] =  - wGp*sgp*DN[i][a];
        lRes[a][i] += - wGp*sgp*DN[i][a];

        // Forcing Term
        u_forcing_term[a][i] += - wGp*fh[i]*lN[iGp][a];
        lRes[a][i] += - wGp*fh[i]*lN[iGp][a];
      }

      // Assemble last d.o.f. for pressure.
      p_grad_term[a] += wGp*trGradHdu*lN[iGp][a]*invEpsilon;
      lRes[a][3]     += wGp*trGradHdu*lN[iGp][a]*invEpsilon;

      p_sg_term[a] += - wGp*sduhDN*invEpsilon;
      lRes[a][3]   += - wGp*sduhDN*invEpsilon;
    }
  }

  // Display residual components for current element
  // Compute the norm for vector quantities
  // u_diff_norm      = femUtils::getMaxModule(u_diff_term, 4, 3);
  // u_conv_norm      = femUtils::getMaxModule(u_conv_term, 4, 3);
  // u_pressgrad_norm = femUtils::getMaxModule(u_pressgrad_term, 4, 3);
  // u_conv_sg_norm   = femUtils::getMaxModule(u_conv_sg_term, 4, 3);
  // u_forcing_norm   = femUtils::getMaxModule(u_forcing_term, 4, 3);
  // p_grad_norm      = femUtils::getMaxModule(p_grad_term);
  // p_sg_norm        = femUtils::getMaxModule(p_sg_term);

  // if(loopTime > 2){
  //   printf("%s %s %s %s %s %s %s\n","u_diff_norm","u_conv_norm","u_pressgrad_norm","u_conv_sg_norm","u_forcing_norm","p_grad_norm","p_sg_norm");
  //   printf("%6.3e %6.3e %6.3e %6.3e %6.3e %6.3e %6.3e\n",u_diff_norm,
  //                                                       u_conv_norm,
  //                                                       u_pressgrad_norm,
  //                                                       u_conv_sg_norm,
  //                                                       u_forcing_norm,
  //                                                       p_grad_norm,
  //                                                       p_sg_norm);
  // }

  // Assemble to global RHS.
  for (uint a = 0; a < 4; ++a){
    for (uint i = 0; i < 4; ++i){
      RHS[nodeIds[a]][i] += lRes[a][i];
    }
  }
}

/* Assemble the projection RHS.
 * elmTaus: (nELms, 1) store each element's tau, update at each time step
 * projRHS: (4*4*NElms, 1) global RHS of projection calculation, needs sync
 * sdus: (4, 3, nElms) velocity subscales evaluated at integration point of each element
 */
void assemble_projRes(const long iElm, 
                      const femModel *model,
                      const femDoubleVec& volumes, 
                      const femDoubleVec& DNs,
                      const femDoubleMat& sol_n, 
                      const femDoubleMat& sol_prev, 
                      const femDoubleVec& char_length, 
                      const femDoubleVec& params, 
                      // return
                      femDoubleVec& elmTaus, 
                      femDoubleVec& sdus, 
                      femDoubleVec& sgps,
                      femDoubleVec& projRHS){

  long nodeIds[4];
  double Ve;
  double DN[3][4]; // 3*4
    
  double sdu[4][3];
  double hdu[4][3]; // 4*3
  double hp[4];

  double du[4][3]; // 4*3
  double p[4];
  double gradDu[3][3];
  double gradP[3];

  // values at Gaussian integration points
  double wGp;
  double hah[3];
  double ahGradDu[3];

  double projL[4][4];
  double proj_den[4][4];
  double maxA = 0.0;
  double tau_u = 0.0;

  // parameters
  double nu = params[1]/params[0];
  double invEps = params[2];
  double c1 = params[3];
  double c2 = params[4];
  
  // Get time step 
  double dt = model->timeStep;

  // inscribe diameter for current element
  double insDiameter = char_length[iElm];

  // Get number of elements
  int nNodes = model->nodeList.size();
  int nElms = model->elementList.size();

  // Divergence of the velocity at the Gauss point
  double divu = 0.0;

  // Remember element's nodeIds.
  for (uint i = 0; i < 4; ++i){
      
    // Get node connections
    nodeIds[i] = model->elementList[iElm]->elementConnections[i];

    // Calculate initial values.
    for (uint j = 0; j < 3; ++j){

      du[i][j] = sol_n[nodeIds[i]][j];
      
      hdu[i][j] = 1.5*du[i][j] - 0.5*sol_prev[nodeIds[i]][j];
          
      // Value of the subgrid velocity at the Gauss points
      sdu[i][j] = sdus[(i*3+j)*nElms+iElm];
      
      projL[i][j] = 0.0;

      proj_den[i][j] = 0.0;

    }      
    p[i] = sol_n[nodeIds[i]][3];
  }

  // Get element volume
  Ve = volumes[iElm];

  // Local element copy of shape function derivatives
  for (uint i = 0; i < 3; ++i){
    for (uint j = 0; j < 4; ++j){
      DN[i][j] = DNs[(i*4+j)*nElms+iElm];
    }
  }

  // Calculate velocity gradient for current element - A CONSTANT FOR TET4
  for (uint i = 0; i < 3; ++i){ // partial u_x, u_y, u_z
    for (uint j = 0; j < 3; ++j){ // partial x, y, z
      gradDu[i][j] = DN[j][0]*du[0][i] \
                   + DN[j][1]*du[1][i] \
                   + DN[j][2]*du[2][i] \
                   + DN[j][3]*du[3][i];
    }
    // Compute pressure gradient for current element - A CONSTANT FOR TET4
    gradP[i] = DN[i][0]*p[0] \
             + DN[i][1]*p[1] \
             + DN[i][2]*p[2] \
             + DN[i][3]*p[3];
  }

  // Gauss point loop
  for (uint iGp = 0; iGp < 4; ++iGp){
    
    // Get GP weights
    wGp = w[iGp] * Ve;

    // Convection at Gaussian point iGp.
    for (uint i = 0; i < 3; ++i){
      // Convective velocity at iGp
      hah[i] = hdu[0][i]*lN[iGp][0] \
             + hdu[1][i]*lN[iGp][1] \
             + hdu[2][i]*lN[iGp][2] \
             + hdu[3][i]*lN[iGp][3] \
             // Sum subgrid velocity already available at the Gauss point
             + sdu[iGp][i];
    }

    // Take max of velocity module over Gauss points
    maxA = fmax(maxA, sqrt(hah[0]*hah[0]+hah[1]*hah[1]+hah[2]*hah[2]));

    // hah dot gradDu - Compute convective term
    divu = 0.0;
    for (uint i = 0; i < 3; ++i){
      ahGradDu[i] = gradDu[i][0]*hah[0] \
                  + gradDu[i][1]*hah[1] \
                  + gradDu[i][2]*hah[2];
      divu += gradDu[i][i];
    }

    // Update subscale at Gauss points for current element - Prepare for subtracting the projected part
    for (uint i = 0; i < 3; ++i){
      sdus[(iGp*3+i)*nElms+iElm] = sdu[iGp][i]/dt - ahGradDu[i] - gradP[i];
    }

    // Update subscale pressures
    sgps[(iGp)*nElms+iElm] += sgps[(iGp)*nElms+iElm]/dt  - invEps*divu;

    // Assemble to local residual.
    for (uint a = 0; a < 4; ++a){
      for (uint i = 0; i < 3; ++i){
        projL[a][i] += wGp*(ahGradDu[i] + gradP[i])*lN[iGp][a];
        // Ve/4 equal at all nodes
        proj_den[a][i] += wGp*lN[iGp][a];
      }
      // Assemble the degree of freedom for the pressure
      projL[a][3] += wGp*(invEps * divu)*lN[iGp][a];
      proj_den[a][3] += wGp*lN[iGp][a];
    }
  }

  // Assemble to global projection RHS.
  // Update the projection for pressure gradient component.
  // NOT SURE I UNDERSTAND - WHY PRESSURE GRADIENT AT THE CENTROID!!!
  for (uint a = 0; a < 4; ++a){    
    for (uint i = 0; i < 4; ++i){
      // Added the Pressure degree of freedom
      // Element-based
      projRHS[(a*4+i)*nElms+iElm] += projL[a][i]/proj_den[a][i];
    }
  }

  // Calculate the elemental tau value.
  tau_u = 1.0/(c1*nu/(insDiameter*insDiameter) + c2*maxA/insDiameter);
  elmTaus[0*nElms + iElm] = 1.0/(1.0/dt + 1.0/tau_u);
  elmTaus[1*nElms + iElm] = (insDiameter*insDiameter)/tau_u;
}

/* 
Update the velocity subscale.
*/
  void update_subscales(bool include_sg_p,
                        const long iElm, 
                        const femModel *model,
                        const femDoubleVec& elmTaus,                        
                        const femDoubleVec& projRHS, 
                        // returns
                        femDoubleVec& sdus,
                        femDoubleVec& sgps){

  long nodeIds[4];
  double projL[4][4];
  double projIGp;
  ulint nElms = model->elementList.size();
  double tau1 = elmTaus[(0)*nElms+iElm];
  double tau2 = elmTaus[(1)*nElms+iElm];
  ulint nNodes = model->nodeList.size();

  // Remember element's nodeIds.
  for (uint i = 0; i < 4; ++i){
    // Get connections
    nodeIds[i] = model->elementList[iElm]->elementConnections[i];

    // Reading projection values at vertices - include pressure
    for (uint j = 0; j < 4; ++j){
      projL[i][j] = projRHS[(i*4+j)*nElms+iElm];
    }
  }

  for (uint iGp = 0; iGp < 4; ++iGp){
    for (uint i = 0; i < 3; ++i){
      projIGp = projL[0][i]*lN[iGp][0] \
              + projL[1][i]*lN[iGp][1] \
              + projL[2][i]*lN[iGp][2] \
              + projL[3][i]*lN[iGp][3];
          
      sdus[(iGp*3+i)*nElms+iElm] = (sdus[(iGp*3+i)*nElms+iElm] + projIGp)*tau1;

    }
    projIGp = projL[0][3]*lN[iGp][0] \
            + projL[1][3]*lN[iGp][1] \
            + projL[2][3]*lN[iGp][2] \
            + projL[3][3]*lN[iGp][3];
    if(include_sg_p){
      sgps[(iGp)*nElms+iElm] = (sgps[(iGp)*nElms+iElm] + projIGp)*tau2;
    }else{
      sgps[(iGp)*nElms+iElm] = 0.0;
    }
    
  }
}

void apply_zeroTraction(const femModel* model,
                        femDoubleMat& sol_n){
  for (ulint i = 0; i < model->presNodesID.size(); i++){
    sol_n[model->presNodesID[i]][3] = model->presNodesValue[i];
  }
}

/* Constant non-reflection outlet B.C. */
void apply_outletBC(const long nNodes, const long nOutlet,
                    const long *outletIndices, const long *outletNghbrNodeIds, 
                    const double *outletNghbrsNs, const double *preDuP, 
                    double *duP){

  // uint idx = get_global_id(0);
  ulint idx = 0;

  const long *iNghbrNodeIds = outletNghbrNodeIds + 4 * idx;

  // Only updates z direction, dof 2. TODO:: use normal!
  double exValue = 0.0; // extrapolation
  for (uint i = 0; i < 4; ++i){
    exValue += outletNghbrsNs[i*nOutlet+idx] * preDuP[2*nNodes+iNghbrNodeIds[i]];
  }
  
  duP[2*nNodes+outletIndices[idx]] = exValue;
}

void apply_linear_outletBC(const long nNodes,
                           const long nOutlet, const long *outletIndices,
                           const long *outletNghbrNodeIds, const double *outletNghbrsNs, 
                           const long *outletNeiNghbrNodeIds, const double *outletNeiNghbrsNs, 
                           const double *preDuP, double *duP){

  // uint idx = get_global_id(0);
  ulint idx = 0;

  const long *iNghbrNodeIds = outletNghbrNodeIds + 4 * idx;
  const long *iSndNghbrNodeIds = outletNeiNghbrNodeIds + 4 * idx;

  // Only updates z direction, dof 2. TODO:: use normal!
  double exValue = 0.0; // extrapolation, first neighbor
  double sndExValue = 0.0; // second neighbor
  for (uint i = 0; i < 4; ++i){
    exValue += outletNghbrsNs[i*nOutlet+idx] * preDuP[2*nNodes+iNghbrNodeIds[i]];
    sndExValue += outletNeiNghbrsNs[i*nOutlet+idx] * preDuP[2*nNodes+iSndNghbrNodeIds[i]];
  }

  // duP[2*nNodes+outletIndices[idx]] = 2.0*exValue - sndExValue;
  duP[2*nNodes+outletIndices[idx]] = 2.0 * sndExValue - exValue;
}

// ===================
// EXPLICIT VMS SOLVER
// ===================
void femGPUFluidSolver::solve(femModel* model){

  // HARDCODED OPTIONS
  bool include_sg_p = true;

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
  femDoubleMat sol_incr;
  femUtils::matZeros(sol_incr,totNodes,model->maxNodeDofs);  
  // ELEMENT-BASED RESIDUAL VECTOR
  femDoubleMat rhs;
  femUtils::matZeros(rhs,totNodes,4);
  // ELEMENT-BASED SUB-GRID SCALE VELOCITY INITIALIZED TO ZERO
  // SUB-GRID VELOCITY IS STORED AT THE 4 GAUSS POINTS
  femDoubleVec sdus(totElements*4*3,0.0);
  // SUB-GRID PRESSURE IS STORED AT THE 4 GAUSS POINTS
  femDoubleVec sgps(totElements*4,0.0);

  // SET INITIAL INFLOW VELOCITIES
  model->setNodeVelocity(currentTime,sol_n);
  // SET ZERO PRESSURES
  apply_zeroTraction(model,sol_n);
  // SET WALL BOUNDARY CONDITIONS
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
  femDoubleVec DNs(totElements*4*3,0.0);
  femDoubleVec char_length(totElements,0.0);
  // Lumped mass matrix - element based
  femDoubleVec nd_mass(totNodes,0.0);
  
  // Assemble initial quantities
  double detJ = 0.0;
  femDoubleMat globShDeriv;
  double curr_vol;
  for(uint loopElement=0;loopElement<totElements;loopElement++){
    // volumes/DNs checked OK
    initial_assemble(loopElement, model,
                     char_length,volumes,DNs,nd_mass);
  }

  // Compute minimum element length for CFL computation
  double min_el_h = std::numeric_limits<double>::max();
  for(uint loopElement=0;loopElement<totElements;loopElement++){
    if(pow(6.0*volumes[loopElement],1.0/3.0) < min_el_h){
      min_el_h = pow(6.0*volumes[loopElement],1.0/3.0);
    }
  }

  // VARIABLES FOR TIME LOOP
  long saveCounter = 0;
  femDoubleMat ha;
  femUtils::matZeros(ha,4,3);
  femDoubleMat sdu;
  femUtils::matZeros(sdu,4,3);
  // Viscosity
  double nu = model->vmsProps[1]/model->vmsProps[0];
  // characteristic length
  double el_h = 0.0;
  // Stabilization constants
  double c1 = model->vmsProps[3];
  double c2 = model->vmsProps[4];
  femDoubleVec elmTaus(totElements*2,0.0);
  // Projected residual - Element-based - Considers velocities and pressures
  femDoubleVec projRHS(totElements*4*4);
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
  double curr_cfl = 0.0;
  double node_mass = 0.0;

  femDoubleVec el_rhs(16,0.0);
  femDoubleMat el_lhs;
  femUtils::matZeros(el_lhs,16,16);
  femDoubleVec sol;
  double curr_sum = 0.0;

  // =========
  // TIME LOOP
  // =========
  for(uint loopTime=1;loopTime<model->totalSteps;loopTime++){

    // Update current time
    currentTime += model->timeStep;

    for(ulint i=0;i<totNodes;i++){
      for(ulint j=0;j<4;j++){
        rhs[i][j] = 0.0;
      }      
    }

    for(ulint loopB=0;loopB<projRHS.size();loopB++){
      projRHS[loopB] = 0.0;
    }
    
    // INITIALIZE ELEMENT RESIDUALS FOR PLOTTING
    for(ulint loopA=0;loopA<totElements;loopA++){
      el_p_res[loopA] = 0.0;
      for(ulint loopB=0;loopB<3;loopB++){
        el_u_res[loopA][loopB] = 0.0;
        el_u_sdus[loopA][loopB] = 0.0;
      }
    }

    // Element Loop
    for(ulint loopElement=0;loopElement<totElements;loopElement++){

      // ASSEMLE ELEMENT RESIDUAL
      assemble_RHS(loopElement,loopTime,model,
                   volumes,DNs,
                   sol_n,sol_prev,
                   // Sub-grid velocities and pressures at all elements and Gauss points
                   sdus,sgps,
                   model->vmsProps,
                   // returns
                   // Local element-based Subgrid velocity at Gauss points
                   sdu,
                   // Local element-based convective velocity at Gauss points
                   ha,
                   // Element-based nodal residual
                   rhs);     

      // COMPUTE RESIDUALS AT ELEMENT CENTROID
      for(ulint loopA=0;loopA<4;loopA++){
        curr_node = model->elementList[loopElement]->elementConnections[loopA];
        for(ulint loopB=0;loopB<3;loopB++){
          el_u_res[loopElement][loopB] += rhs[curr_node][loopB]*0.25;
        }
        el_p_res[loopElement] += rhs[curr_node][3]*0.25;
      }

      // COMPUTE SUBGRID VELOCITY AT ELEMENT CENTROID - SDU computed at Gauss Points!!!
      for(ulint loopA=0;loopA<4;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          el_u_sdus[loopElement][loopB] += sdu[loopA][loopB]*0.25;
        }
      }
      
      // Compute Projected Residual
      assemble_projRes(loopElement,model,
                       volumes,DNs,
                       sol_n,sol_prev, 
                       char_length, model->vmsProps,
                       // return
                       // Element stabilization coefficients
                       elmTaus, 
                       // Updated Subgrid scale with residual subtracted
                       sdus, 
                       // Updated Subgrid pressure with residual subtracted
                       sgps,
                       // Projected residual
                       projRHS);

      // Update Subgrid scale velocity
      update_subscales(include_sg_p,
                       loopElement, 
                       model,
                       elmTaus,
                       projRHS, 
                       // returns
                       sdus,
                       sgps);

    } // End of Element loop

    for(ulint i=0;i<totNodes;i++){
      for(ulint j=0;j<4;j++){
        sol_n[i][j]  = sol_prev[i][j] - model->timeStep * rhs[i][j] / nd_mass[i];
      }      
    }

    // COMPUTE MAX RESIDUAL NORMS
    el_u_res_norm = 0.0;
    el_p_res_norm = 0.0;
    for(ulint loopElement=0;loopElement<totElements;loopElement++){
      curr_u_res = sqrt(el_u_res[loopElement][0]*el_u_res[loopElement][0] 
                      + el_u_res[loopElement][1]*el_u_res[loopElement][1] 
                      + el_u_res[loopElement][2]*el_u_res[loopElement][2]);
      if(curr_u_res > el_u_res_norm){
        el_u_res_norm = curr_u_res;
      }
      if(fabs(el_p_res[loopElement]) > el_p_res_norm){
        el_p_res_norm = fabs(el_p_res[loopElement]);
      }
    }

    // APPLY INLET VELOCITIES
    model->setNodeVelocity(currentTime,sol_n);
    // APPLY OUTLET BC
    apply_zeroTraction(model,sol_n);
    
    // APPLY OUTLET BC - NO OUTLET BC FOR NOW
    // apply_outletBC(nNodes,nOutlet,
    //                 outletIndices, outletNghbrNodeIds, 
    //                 const double *outletNghbrsNs, const double *preDuP, 
    //                 double *duP);

    // apply_linear_outletBC(nNodes,
    //                       nOutlet, outletIndices,
    //                       outletNghbrNodeIds, outletNghbrsNs, 
    //                       outletNeiNghbrNodeIds, outletNeiNghbrsNs, 
    //                       preDuP, duP);

    // SET DIRICHLET BC
    model->setDirichletBC(sol_n);

    // GET CFL
    // WRITE CFL MESSAGE
    // Get maximum velocity 
    double max_ini_el_u = 0.0;
    for(uint loopNode=0;loopNode<totNodes;loopNode++){
      if(sqrt(sol_n[loopNode][0]*sol_n[loopNode][0] + sol_n[loopNode][1]*sol_n[loopNode][1] + sol_n[loopNode][2]*sol_n[loopNode][2]) > max_ini_el_u){
        max_ini_el_u = sqrt(sol_n[loopNode][0]*sol_n[loopNode][0] + sol_n[loopNode][1]*sol_n[loopNode][1] + sol_n[loopNode][2]*sol_n[loopNode][2]);
      }
    }
    curr_cfl = max_ini_el_u*model->timeStep/min_el_h;

    // WRITE MESSAGES
    if(loopTime == 0){
      printf("--- STARTING TIME LOOP\n");
      printf("%9s %20s %20s %20s %10s\n","TIME LOOP","TIME [s]","MAX RES NORM U","MAX RES NORM P","CFL");
    }
    printf("%9d %20.8f %20.3e %20.3e %10.3f\n",loopTime,currentTime,el_u_res_norm,el_p_res_norm,curr_cfl);

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



