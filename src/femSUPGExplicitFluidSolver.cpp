#include "femSUPGExplicitFluidSolver.h"
                          
// CONSTRUCTOR
femSUPGExplicitFluidSolver::femSUPGExplicitFluidSolver(){
}

void write_log_header(string file_name){
    // Append to file
    FILE* fptr = fopen(file_name.c_str(), "a");
    // Write some text to the file
    fprintf(fptr, "%5s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n","IT",
                                                                                          "Time",
                                                                                          "min X","max X",
                                                                                          "min Y","max Y",
                                                                                          "min Z","max Z",
                                                                                          "min tau","max tau",
                                                                                          "min lambda","max lambda",
                                                                                          "CFL v",
                                                                                          "CFL a",
                                                                                          "E_min","E_max",
                                                                                          "nu_min","nu_max");
    // Close the file
    fclose(fptr); 
    // Write to screen
    printf("%5s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n","IT",
                                                                                          "Time",
                                                                                          "min X","max X",
                                                                                          "min Y","max Y",
                                                                                          "min Z","max Z",
                                                                                          "min tau","max tau",
                                                                                          "min lamda","max lamda",
                                                                                          "CFL v",
                                                                                          "CFL a",
                                                                                          "E_min","E_max",
                                                                                          "nu_min","nu_max");
}


// ASSEMBLE IN GLOBAL STIFFNESS MATRIX
void assemble_in_glob(femModel* model, const int loopElement, const femDoubleMat& K_mat,femDoubleMat& glob_K){
  ulint curr_node_row = 0;
  ulint curr_node_col = 0;
  for(int loopA=0;loopA<model->elementList[loopElement]->elementConnections.size();loopA++){
    curr_node_row = model->elementList[loopElement]->elementConnections[loopA];
    for(int loopB=0;loopB<model->elementList[loopElement]->elementConnections.size();loopB++){
      curr_node_col = model->elementList[loopElement]->elementConnections[loopB];
      for(int loopC=0;loopC<3;loopC++){
        for(int loopD=0;loopD<3;loopD++){
          glob_K[curr_node_row*3 + loopC][curr_node_col*3 + loopD] += K_mat[loopA*3 + loopC][loopB*3 + loopD];
        }
      }
    }
  }
}

// APPLY DIRICHLET CONDITIONS TO MATRIX
void apply_dir_vel_glob_K_f(femModel* model,const double currTime,femDoubleMat& glob_K,femDoubleVec& glob_f){

  int curr_node = 0;
  femDoubleVec curr_vel;
  femDoubleMat table;

  // APPLY DIRICHLET VALUES
  for(ulint loopA=0;loopA<model->diricheletBCNode.size();loopA++){
    curr_node = model->diricheletBCNode[loopA];
    for(ulint loopB=0;loopB<3;loopB++){
      // SET ROWS TO UNIT DIAGONAL
      for(ulint loopC=0;loopC<glob_K[curr_node*3 + loopB].size();loopC++){
        if(loopC == curr_node*3 + loopB){
          glob_K[curr_node*3 + loopB][loopC] = 1.0;
        }else{
          glob_K[curr_node*3 + loopB][loopC] = 0.0;
        }
      }        
      // SET RHS TO DIRICHLET VALUES
      glob_f[curr_node*3 + loopB] = model->diricheletBCValues[loopA][loopB];
    }
  }

  // Loop over the nodes with prescribed velocities
  for(ulint loopA=0;loopA<model->velNodesID.size();loopA++){ 
    curr_node = model->velNodesID[loopA];        
    table = model->velNodesTimeVals[loopA];
    curr_vel = femUtils::interpTableData(currTime, table);
    for(ulint loopB=0;loopB<3;loopB++){
      for(ulint loopC=0;loopC<glob_K[curr_node*3 + loopB].size();loopC++){
        if(loopC == curr_node*3 + loopB){
          glob_K[curr_node*3 + loopB][loopC] = 1.0;
        }else{
          glob_K[curr_node*3 + loopB][loopC] = 0.0;
        }
      }        
      // SET RHS TO DIRICHLET VALUES
      glob_f[curr_node*3 + loopB] = curr_vel[loopB];
    }
  }
}

void save_inital_solution(femModel* model, 
                          const femDoubleMat& sol,
                          const femDoubleVec& pres){

  // CREATE MODEL RESULTS
  // ADD VELOCITY
  femResult* res = new femResult();
  femDoubleVec res_temp;
  res->label = string("velocity");
  res->type = frNode;
  res->numComponents = 3;
  // Assign to values
  for(size_t loopA=0;loopA<sol.size();loopA++){
    res_temp.clear();
    res_temp.push_back(sol[loopA][0]);
    res_temp.push_back(sol[loopA][1]);
    res_temp.push_back(sol[loopA][2]);
    res->values.push_back(res_temp);
  }
  model->resultList.push_back(res);
  // ADD PRESSURE
  res = new femResult();
  res->label = string("pressure");
  res->type = frElement;
  res->numComponents = 1;
  // Assign to values
  for(size_t loopA=0;loopA<pres.size();loopA++){
    res_temp.clear();
    res_temp.push_back(pres[loopA]);
    res->values.push_back(res_temp);
  }
  model->resultList.push_back(res);
  model->ExportToVTKLegacy(string(model->sol_file_prefix + "_" + femUtils::intToStr(0) + ".vtk"),false);
}

// =============
// EVAL L MATRIX
// =============
void eval_L_matrix(const femDoubleVec& shapeFunction,
                    femDoubleMat& L_mat){
                    
  // Make sure the matrix is initialized
  for(ulint loopA=0;loopA<L_mat.size();loopA++){                    
    for(ulint loopB=0;loopB<L_mat[loopA].size();loopB++){
      L_mat[loopA][loopB] = 0.0;
    }
  }
  // Fill matrix
  for(ulint loopA=0;loopA<shapeFunction.size();loopA++){
    L_mat[0][loopA*3 + 0] = shapeFunction[loopA];
    L_mat[1][loopA*3 + 1] = shapeFunction[loopA];
    L_mat[2][loopA*3 + 2] = shapeFunction[loopA];
  }
}

// =============
// EVAL U MATRIX
// =============
void eval_U_matrix(const femDoubleVec& u_gp,
                   femDoubleMat& U_mat){

  // Make sure the matrix is initialized
  for(ulint loopA=0;loopA<U_mat.size();loopA++){
    for(ulint loopB=0;loopB<U_mat[loopA].size();loopB++){
      U_mat[loopA][loopB] = 0.0;
    }
  }
  // Compute the velocity compoments at the current Gauss point
  // First row
  U_mat[0][0] = u_gp[0];
  U_mat[0][1] = u_gp[1];
  U_mat[0][2] = u_gp[2];
  // Second row
  U_mat[1][3] = u_gp[0];
  U_mat[1][4] = u_gp[1];
  U_mat[1][5] = u_gp[2];
  // Third row
  U_mat[2][6] = u_gp[0];
  U_mat[2][7] = u_gp[1];
  U_mat[2][8] = u_gp[2];
}


// =============
// EVAL S MATRIX
// =============
void eval_S_matrix(const femDoubleMat& shapeDeriv,
                   femDoubleMat& S_mat){

  for(ulint loopA=0;loopA<S_mat.size();loopA++){
    for(ulint loopB=0;loopB<S_mat[loopA].size();loopB++){
      S_mat[loopA][loopB] = 0.0;
    }
  }

  for(ulint loopA=0;loopA<shapeDeriv.size();loopA++){

    S_mat[0][loopA*3 + 0] = shapeDeriv[loopA][0];
    S_mat[1][loopA*3 + 0] = shapeDeriv[loopA][1];
    S_mat[2][loopA*3 + 0] = shapeDeriv[loopA][2];

    S_mat[3][loopA*3 + 1] = shapeDeriv[loopA][0];
    S_mat[4][loopA*3 + 1] = shapeDeriv[loopA][1];
    S_mat[5][loopA*3 + 1] = shapeDeriv[loopA][2];

    S_mat[6][loopA*3 + 2] = shapeDeriv[loopA][0];
    S_mat[7][loopA*3 + 2] = shapeDeriv[loopA][1];
    S_mat[8][loopA*3 + 2] = shapeDeriv[loopA][2];

  }
}

// Eval B matrix
void eval_B_matrix(const bool use_b_hat,
                   const femDoubleMat& shapeDeriv,
                   // Derivative of the shape functions at element centroid
                   const femDoubleMat& shapeDeriv_0, 
                   femDoubleMat& B_mat){

  for(ulint loopA=0;loopA<B_mat.size();loopA++){
    for(ulint loopB=0;loopB<B_mat[loopA].size();loopB++){
      B_mat[loopA][loopB] = 0.0;
    }
  }

  double B_1 = 0.0;
  double B_2 = 0.0;
  double B_3 = 0.0;
  double B_4 = 0.0;
  double B_5 = 0.0;
  double B_6 = 0.0;

  for(ulint loopA=0;loopA<shapeDeriv.size();loopA++){

    if(use_b_hat){
      // INCOMPRESSIBLE FORMULATION
      // Compute the modified volumetric terms
      B_1 = (1.0/3.0) * (shapeDeriv_0[loopA][0] - shapeDeriv[loopA][0]);
      B_2 = shapeDeriv[loopA][0] - B_1;
      B_3 = (1.0/3.0) * (shapeDeriv_0[loopA][1] - shapeDeriv[loopA][1]);
      B_4 = shapeDeriv[loopA][1] - B_3;
      B_5 = (1.0/3.0) * (shapeDeriv_0[loopA][2] - shapeDeriv[loopA][2]);
      B_6 = shapeDeriv[loopA][2] - B_5;
      // Assemble nodal contribution
      // Volumetric
      // First row
      B_mat[0][loopA*3 + 0] = B_2;
      B_mat[0][loopA*3 + 1] = B_3;
      B_mat[0][loopA*3 + 2] = B_5;
      // Second row
      B_mat[1][loopA*3 + 0] = B_1;
      B_mat[1][loopA*3 + 1] = B_4;
      B_mat[1][loopA*3 + 2] = B_5;
      // Third row
      B_mat[2][loopA*3 + 0] = B_1;
      B_mat[2][loopA*3 + 1] = B_3;
      B_mat[2][loopA*3 + 2] = B_6;
    }else{
      // COMPRESSIBLE FORMULATION
      // First row
      B_mat[0][loopA*3 + 0] = shapeDeriv[loopA][0];
      B_mat[0][loopA*3 + 1] = 0.0;
      B_mat[0][loopA*3 + 2] = 0.0;
      // Second row
      B_mat[1][loopA*3 + 0] = 0.0;
      B_mat[1][loopA*3 + 1] = shapeDeriv[loopA][1];
      B_mat[1][loopA*3 + 2] = 0.0;
      // Third row
      B_mat[2][loopA*3 + 0] = 0.0;
      B_mat[2][loopA*3 + 1] = 0.0;
      B_mat[2][loopA*3 + 2] = shapeDeriv[loopA][2];
    }
    // Deviatoric
    // Fourth row
    B_mat[3][loopA*3 + 0] = shapeDeriv[loopA][1];
    B_mat[3][loopA*3 + 1] = shapeDeriv[loopA][0];
    B_mat[3][loopA*3 + 2] = 0.0;
    // Fifth row
    B_mat[4][loopA*3 + 0] = 0.0;
    B_mat[4][loopA*3 + 1] = shapeDeriv[loopA][2];
    B_mat[4][loopA*3 + 2] = shapeDeriv[loopA][1];
    // Sixth row
    B_mat[5][loopA*3 + 0] = shapeDeriv[loopA][2];
    B_mat[5][loopA*3 + 1] = 0.0;
    B_mat[5][loopA*3 + 2] = shapeDeriv[loopA][0];
  }
}

// =============
// EVAL D MATRIX
// =============
void eval_D_matrix(double lam,double mu,femDoubleMat& D_mat){

  for(ulint loopA=0;loopA<D_mat.size();loopA++){
    for(ulint loopB=0;loopB<D_mat[loopA].size();loopB++){
      D_mat[loopA][loopB] = 0.0;
    }
  }

  D_mat[0][0] = lam + 2*mu;
  D_mat[0][1] = lam;
  D_mat[0][2] = lam;

  D_mat[1][0] = lam;
  D_mat[1][1] = lam + 2*mu;
  D_mat[1][2] = lam;

  D_mat[2][0] = lam;
  D_mat[2][1] = lam;
  D_mat[2][2] = lam + 2*mu;

  D_mat[3][3] = mu;
  D_mat[4][4] = mu;
  D_mat[5][5] = mu;
}

// =============
// EVAL H MATRIX
// =============
void eval_H_matrix(const femDoubleMat& shapeDeriv,
                   femDoubleVec& H_mat){

  for(ulint loopA=0;loopA<H_mat.size();loopA++){
      H_mat[loopA] = 0.0;
  }

  for(ulint loopA=0;loopA<shapeDeriv.size();loopA++){
    H_mat[loopA*3 + 0] = shapeDeriv[loopA][0];
    H_mat[loopA*3 + 1] = shapeDeriv[loopA][1];
    H_mat[loopA*3 + 2] = shapeDeriv[loopA][2];
  }
}

// ================================================
// COMPUTE STABILIZATION COEFFICIENT AT GAUSS POINT
// ================================================
void eval_stabilization(const femDoubleVec& u_gp,
                        double rho, double mu, double cI, double dt,
                        const femDoubleMat& elGeomMat,
                        double& tau_SUPG,
                        double& nu_C){

  double tau_term2 = 0.0;
  double tau_term3 = 0.0;
  double temp = 0.0;
  double trG = elGeomMat[0][0] + elGeomMat[1][1] + elGeomMat[2][2];

  // Timestep dependent stabilization term
  double tau_term1 = (2.0/dt)*(2.0/dt);

  // Velocity-dependent stabilization term
  tau_term2 = 0.0;
  for(ulint loopA=0;loopA<3;loopA++){
    temp = 0.0;
    for(ulint loopB=0;loopB<3;loopB++){
      temp += elGeomMat[loopA][loopB] * u_gp[loopB];
    }
    tau_term2 += u_gp[loopA] * temp;
  }

  // Third term of tau SUPG  
  tau_term3 = 0.0;
  for(ulint loopA=0;loopA<3;loopA++){
    for(ulint loopB=0;loopB<3;loopB++){
      tau_term3 += elGeomMat[loopA][loopB] * elGeomMat[loopA][loopB];
    }
  }
  tau_term3 = tau_term3 * (mu/rho) * (mu/rho) * cI;

  // printf("Stabilization terms: 1 %f, 2 %f, 3 %f\n",tau_term1,tau_term2,tau_term3);

  // Assemble tau_SUPG
  // tau_SUPG = 1.0/sqrt(tau_term1 + tau_term2 + tau_term3);
  // CAREFUL!!! RESTORE TIME STEP TERM!!!
  tau_SUPG = 1.0/sqrt(tau_term2 + tau_term3);
  // Assemble n_C
  nu_C = 1.0/(trG * tau_SUPG);  
}

// ===================================
// COMPUTE BULK MODULUS AT GAUSS POINT
// ===================================
void eval_bulk_modulus(const femDoubleVec& u_gp,
                        double rho, double mu, double alpha, 
                        double dt,
                        const femDoubleMat& elGeomMat,
                        double l_X,double l_Y,double l_Z,
                        double& lam){

    femDoubleMat invG;
    femUtils::matZeros(invG,3,3);
    double detJ = 0.0;
    femDoubleVec invGu(3,0.0);

    // Compute inverse of geometric matrix
    femUtils::invert3x3Matrix(elGeomMat,invG,detJ);

    // Compute the velocity norm
    double u_norm = femUtils::DoEucNorm(u_gp);

    // Determine the first limit
    double limit_1 = (mu/rho);

    // Determine the second limit
    for(int loopA=0;loopA<3;loopA++){
      invGu[loopA] = 0.0;
      for(int loopB=0;loopB<3;loopB++){
        invGu[loopA] += invG[loopA][loopB] * u_gp[loopB];
      }
    }
    double limit_2 = 0.0;
    for(int loopA=0;loopA<3;loopA++){
      limit_2 += invGu[loopA] * u_gp[loopA];
    }
    limit_2 = sqrt(limit_2);

    // Compare second limit
    // printf("velocity: %e %e %e\n",u_gp[0],u_gp[1],u_gp[2]);
    // printf("lam: second limit: %e %e %e %e\n",limit_2,u_norm*l_X/2.0,u_norm*l_Y/2.0,u_norm*l_Z/2.0);

    // Determine the third limit
    double limit_3 = 0.0;
    double trace_inG = 0.0;
    for(int loopA=0;loopA<3;loopA++){
      trace_inG += invG[loopA][loopA];
    }
    limit_3 = (trace_inG)/(alpha * dt);
    // printf("lam: third limit: %e %e %e %e\n",trace_inG,l_X*l_X,l_Y*l_Y,l_Z*l_Z);
    // getchar();

    // Compute the resulting alpha
    // lam = alpha * rho * max(limit_1,max(limit_2,limit_3));
    lam = alpha * rho * max(limit_1,limit_2);
}

// ===============================
// COMPUTE POST_PROCESSED PRESSURE
// ===============================
void compute_pressure(femModel* model, const femDoubleVec& el_lam,const femDoubleMat& sol_n, femDoubleVec& pres_n){

  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);
  femDoubleVec el_sol(8*3,0.0);
  double detJ = 0.0;
  femDoubleMat shapeDeriv;
  femDoubleVec H_vec(8*3,0.0);

  for(ulint loopElement=0;loopElement<model->elementList.size();loopElement++){

      // Get total number of nodes in element
      femIntVec el_conn = model->elementList[loopElement]->elementConnections;
      ulint tot_el_nodes = el_conn.size();

      // Get solution vector for current element
      for(ulint loopA=0;loopA<tot_el_nodes;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          el_sol[loopA*3 + loopB] = sol_n[el_conn[loopA]][loopB];
        }
      }

      // Get Gauss Points and Weights
      femDoubleMat gps = rule->getCoords(tot_el_nodes,d3);
      femDoubleVec gpw = rule->getWeights(tot_el_nodes,d3);

      // GAUSS POINT LOOP
      double numerator = 0.0;
      double denominator = 0.0;
      double div_u_igp = 0.0;
      for(int igp=0;igp<rule->getTotGP(tot_el_nodes,d3);igp++){

        // Eval Current Shape Derivatives Matrix
        model->elementList[loopElement]->evalGlobalShapeFunctionDerivative(model->nodeList,
                                                                           gps[igp][0],gps[igp][1],gps[igp][2],
                                                                           detJ,shapeDeriv);

        // Eval H matrix
        eval_H_matrix(shapeDeriv,H_vec);

        // Compute mean divergence of the velocity field
        div_u_igp = 0.0;
        for(int loopA=0;loopA<tot_el_nodes*3;loopA++){
          div_u_igp += H_vec[loopA] * el_sol[loopA];
        }

        // integrate the divergence over the element
        numerator += div_u_igp * detJ * gpw[igp];
        denominator += detJ * gpw[igp];

      } // END OF GAUSS POINT LOOP

      // Assign pressure
      pres_n[loopElement] = -el_lam[loopElement] * (numerator/denominator);
      
  } // END OF ELEMENT LOOP  
}

// =========================================
// PRELIMIARY ASSEMBLY OF LUMPED MASS MATRIX
// =========================================

void form_lumped_mass(femModel* model, double rho, femDoubleVec& l_M){

  // fem Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);
  femDoubleVec shapeFunction;
  femDoubleMat shapeDeriv;
  double detJ = 0.0;
  double temp_1 = 0.0;
  femDoubleMat L_mat;
  femUtils::matZeros(L_mat,3,8*3);
  femDoubleVec el_lumped_mat(8*3,0.0);
  femDoubleMat BMass_mat;
  femUtils::matZeros(BMass_mat,8*3,8*3);
  
  // Init mass
  for(ulint loopA=0;loopA<l_M.size();loopA++){
    l_M[loopA] = 0.0;
  }
  
  for(ulint loopElement=0;loopElement<model->elementList.size();loopElement++){

      // Get total number of nodes in element
      femIntVec el_conn = model->elementList[loopElement]->elementConnections;
      ulint tot_el_nodes = el_conn.size();

      // Get Gauss Points and Weights
      femDoubleMat gps = rule->getCoords(tot_el_nodes,d3);
      femDoubleVec gpw = rule->getWeights(tot_el_nodes,d3);

      for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
        for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
          BMass_mat[loopA][loopB] = 0.0;
        }
      }

      // GAUSS POINT LOOP
      for(int igp=0;igp<rule->getTotGP(tot_el_nodes,d3);igp++){


        // Eval Shape Function
        model->elementList[loopElement]->evalShapeFunction(model->nodeList,
                                                           gps[igp][0],gps[igp][1],gps[igp][2],
                                                           shapeFunction);

        // Eval Current Shape Derivatives Matrix
        model->elementList[loopElement]->evalGlobalShapeFunctionDerivative(model->nodeList,
                                                                           gps[igp][0],gps[igp][1],gps[igp][2],
                                                                           detJ,shapeDeriv);

        // Eval L matrix
        eval_L_matrix(shapeFunction,L_mat);

        // FOR MASS, STIFFNESS AND STABILIZATION MATRICES
        for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
          for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
            
            // Bubnov Mass Matrix
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<3;loopC++){
              temp_1 += rho * L_mat[loopC][loopA] * L_mat[loopC][loopB] * detJ * gpw[igp];
            }
            BMass_mat[loopA][loopB] += temp_1;
            
          }
        } // END OF MATRIX COMPONENT LOOP

      } // END OF GAUSS POINT LOOP
      
      //printf("Mass matrix norm: %f\n",femUtils::getMatrixNorm(BMass_mat));
      //getchar();

      // LUMP BUBNOV MASS MATRIX
      for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
        temp_1 = 0.0;
        for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
          temp_1 += BMass_mat[loopA][loopB];
        }
        el_lumped_mat[loopA] = temp_1;
      }

      // printf("%f\n",rho);
      // femUtils::printVector(el_lumped_mat);
      // getchar();

      // ASSEMBLE IN GLOBAL LUMPED MATRIX
      for(ulint loopA=0;loopA<tot_el_nodes;loopA++){
        for(ulint loopB=0;loopB<3;loopB++){
          l_M[el_conn[loopA]*3 + loopB] += el_lumped_mat[loopA*3 + loopB];
        }
      }
  } // END OF ELEMENT LOOP  
}


// ==========================
// GET MAX AND MIN VELOCITIES
// ==========================
void get_solution_extrema(const femDoubleMat& sol_n,
                          double& min_vel_mod,double& max_vel_mod,
                          double& min_vel_X,double& max_vel_X,
                          double& min_vel_Y,double& max_vel_Y,
                          double& min_vel_Z,double& max_vel_Z){
  // Find min and max velocities
  max_vel_mod = 0.0;
  min_vel_mod = std::numeric_limits<double>::max();

  max_vel_X = -std::numeric_limits<double>::max();
  min_vel_X = std::numeric_limits<double>::max();

  max_vel_Y = -std::numeric_limits<double>::max();
  min_vel_Y = std::numeric_limits<double>::max();

  max_vel_Z = -std::numeric_limits<double>::max();
  min_vel_Z = std::numeric_limits<double>::max();

  double curr_mod = 0.0;

  for(ulint loopA=0;loopA<sol_n.size();loopA++){
    curr_mod = sqrt(sol_n[loopA][0]*sol_n[loopA][0] + sol_n[loopA][1]*sol_n[loopA][1] + sol_n[loopA][2]*sol_n[loopA][2]);
    // Module
    if(curr_mod > max_vel_mod){
      max_vel_mod = curr_mod;
    }
    if(curr_mod < min_vel_mod){
      min_vel_mod = curr_mod;
    }
    // X
    if(sol_n[loopA][0] > max_vel_X){
      max_vel_X = sol_n[loopA][0];
    }
    if(sol_n[loopA][0] < min_vel_X){
      min_vel_X = sol_n[loopA][0];
    }
    // Y
    if(sol_n[loopA][1] > max_vel_Y){
      max_vel_Y = sol_n[loopA][1];
    }
    if(sol_n[loopA][1] < min_vel_Y){
      min_vel_Y = sol_n[loopA][1];
    }
    // Z
    if(sol_n[loopA][2] > max_vel_Z){
      max_vel_Z = sol_n[loopA][2];
    }
    if(sol_n[loopA][2] < min_vel_Z){
      min_vel_Z = sol_n[loopA][2];
    }
  }
}

// =====================================================
// EVAL AVERAGE SHAPE DERIVATIVE MATRIX OVER ONE ELEMENT
// =====================================================
void eval_avg_shapeDeriv(femModel* model, ulint loopElement, femDoubleMat& avg_sd){
  
  // fem Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);
  femDoubleMat shapeDeriv;
  femDoubleMat den;
  femUtils::matZeros(den,8,3);

  // Get number of nodes in current element
  int tot_el_nodes = model->elementList[loopElement]->elementConnections.size();

  // Get Gauss Points and Weights
  femDoubleMat gps = rule->getCoords(tot_el_nodes,d3);
  femDoubleVec gpw = rule->getWeights(tot_el_nodes,d3);
  double detJ = 0.0;

  // INIT AVG_SD and DEN
  for(int loopA=0;loopA<8;loopA++){
    for(int loopB=0;loopB<3;loopB++){
      avg_sd[loopA][loopB] = 0.0;
      den[loopA][loopB] = 0.0;
    }
  }

  // GAUSS POINT LOOP
  // printf("Total Gauss Points: %d\n",rule->getTotGP(tot_el_nodes,d3));  
  for(int igp=0;igp<rule->getTotGP(tot_el_nodes,d3);igp++){

    // Eval Current Shape Derivatives Matrix
    model->elementList[loopElement]->evalGlobalShapeFunctionDerivative(model->nodeList,
                                                                        gps[igp][0],gps[igp][1],gps[igp][2],
                                                                        detJ,shapeDeriv);
    // Sum contributions for all components
    for(ulint loopA=0;loopA<8;loopA++){
      for(ulint loopB=0;loopB<3;loopB++){
        avg_sd[loopA][loopB] += shapeDeriv[loopA][loopB] * detJ * gpw[igp];
        den[loopA][loopB] += detJ * gpw[igp];
      }
    }
  }

  // Sum contributions for all components
  for(ulint loopA=0;loopA<8;loopA++){
    for(ulint loopB=0;loopB<3;loopB++){
      avg_sd[loopA][loopB] = avg_sd[loopA][loopB]/den[loopA][loopB];
    }
  }
}

// ==================================
// GET ELEMENT CHARACTERISTIC LENGTHS
// ==================================
void get_el_lengths(femModel* model,int curr_el,double& len_X,double& len_Y,double& len_Z){
  double min_X = std::numeric_limits<double>::max(); double max_X = -std::numeric_limits<double>::max();
  double min_Y = std::numeric_limits<double>::max(); double max_Y = -std::numeric_limits<double>::max();
  double min_Z = std::numeric_limits<double>::max(); double max_Z = -std::numeric_limits<double>::max();
  ulint curr_node = 0;
  double curr_X = 0.0; double curr_Y = 0.0; double curr_Z = 0.0; 
  // Get max and min coords in all directions
  for(int loopA=0;loopA<model->elementList[curr_el]->elementConnections.size();loopA++){
    curr_node  = model->elementList[curr_el]->elementConnections[loopA];
    curr_X = model->nodeList[curr_node]->coords[0];
    curr_Y = model->nodeList[curr_node]->coords[1];
    curr_Z = model->nodeList[curr_node]->coords[2];
    // X Coords
    if(curr_X < min_X){
      min_X = curr_X;
    }
    if(curr_X > max_X){
      max_X = curr_X;
    }
    // Y Coords
    if(curr_Y < min_Y){
      min_Y = curr_Y;
    }
    if(curr_Y > max_Y){
      max_Y = curr_Y;
    }
    // Z Coords
    if(curr_Z < min_Z){
      min_Z = curr_Z;
    }
    if(curr_Z > max_Z){
      max_Z = curr_Z;
    }
  }
  // Compute the lengths
  len_X = max_X - min_X;
  len_Y = max_Y - min_Y;
  len_Z = max_Z - min_Z;
}

// ===================
// EXPLICIT VMS SOLVER
// ===================
void femSUPGExplicitFluidSolver::solve(femModel* model){

  // GET SOLVER OPTIONS
  bool use_euler  = model->exOpts[0];
  bool use_B_hat  = model->exOpts[1];
  bool include_K1 = model->exOpts[2];
  bool include_K2 = model->exOpts[3];
  bool include_K3 = model->exOpts[4];
  bool include_K4 = model->exOpts[5];

  // USE CONSTANT LAM EQUAL TO ALPHA
  bool use_constant_lamda = false;

  // GET LIMIT ON 
  double max_allowed_vel = model->exit_condition;

  // GET GLOBAL FIXED DOFs
  string global_fixed_dof = model->glob_fixed_dofs;

  // OPTIONS
  // TYPE OF EXPLICIT INTEGRATION
  explicit_alg time_alg;
  if(use_euler){
    time_alg = eaForwardEuler;
  }else{
    time_alg = eaAdamsBashforth;
  }

  // SOLVE STATIC PROBLEM FOR DEBUG
  bool solve_static_K = false;

  // SET UP INITIAL CONDITIONS
  double currentTime = 0.0;
  
  // TOTAL QTYs
  long totNodes = model->nodeList.size();
  long totElements = model->elementList.size();
  int max_dofs = model->maxNodeDofs;

  // Get fluid properties
  double rho = model->vmsProps[0];
  double mu  = model->vmsProps[1];
  double alpha = model->vmsProps[2];
  double cI  = model->vmsProps[3];

  // Time step counter
  ulint saveCounter = 0;

  // fem Integration Rule
  femIntegrationRule* rule = new femIntegrationRule(irSecondOrder);

  // Results
  femResult* res;
  femDoubleVec res_temp;

  // AUX
  ulint curr_node = 0.0;
  double temp_1 = 0.0;
  double temp_2 = 0.0;

  // Declare main variables - Set to zero
  // Main solution at t_n+1
  femDoubleMat sol_n;
  femUtils::matZeros(sol_n,totNodes,model->maxNodeDofs);  
  // Post-processed Pressure at t_n+1
  femDoubleVec pres_n(totElements,0.0);
  // Previous soluton at t_n
  femDoubleMat sol_p;
  femUtils::matZeros(sol_p,totNodes,model->maxNodeDofs);
  // Global lumped mass
  femDoubleVec l_M(totNodes*3,0.0);

  // TEST ASSEMBLE GLOBAL STIFFNESS
  femDoubleMat glob_K;
  femDoubleVec glob_f;
  femDoubleVec lin_sol; 
  if(solve_static_K){
    femUtils::matZeros(glob_K,totNodes*model->maxNodeDofs,totNodes*model->maxNodeDofs);
    femDoubleVec glob_f(totNodes*model->maxNodeDofs,0.0);
  }

  // FORCE VECTORS FOR TIME INTEGRATION
  // time n+1
  femDoubleMat k_sol_p;
  femUtils::matZeros(k_sol_p,8*3,totElements);  
  // time n
  femDoubleMat k_sol_pp;
  femUtils::matZeros(k_sol_pp,8*3,totElements);

  // MEAN VALUE OF THE BULK MODULUS OVER THE ELEMENTS
  femDoubleVec el_lam(totElements,0.0);

  // SET WALL BOUNDARY CONDITIONS
  model->setDirichletBC(sol_p);
  // SET INITIAL INFLOW VELOCITIES
  model->setNodeVelocity(currentTime,sol_p);  

  // POST-PROCESS PRESSURE
  compute_pressure(model,el_lam,sol_p,pres_n);

  // SAVE INITIAL CONDITIONS
  save_inital_solution(model,sol_p,pres_n);

  // ASSEMBLE ELEMENT-BASED QTY
  femDoubleMat L_mat;
  femUtils::matZeros(L_mat,3,8*3);  
  
  femDoubleMat U_mat;
  femUtils::matZeros(U_mat,3,9);

  femDoubleMat S_mat;
  femUtils::matZeros(S_mat,9,8*3);

  femDoubleMat US_mat;
  femUtils::matZeros(US_mat,3,8*3);

  femDoubleMat B_mat;
  femUtils::matZeros(B_mat,6,8*3);

  femDoubleMat D_mat;
  femUtils::matZeros(D_mat,6,6);

  femDoubleMat DB_mat;
  femUtils::matZeros(DB_mat,6,8*3);
  
  // INIT LOCAL MATRICES
  femDoubleMat Mass_mat;
  femUtils::matZeros(Mass_mat,8*3,8*3);
  femDoubleMat K_mat;
  femUtils::matZeros(K_mat,8*3,8*3);
  femDoubleMat C_mat;
  femUtils::matZeros(C_mat,8*3,8*3);

  femDoubleMat K1_mat;
  femUtils::matZeros(K1_mat,8*3,8*3);
  femDoubleMat K2_mat;
  femUtils::matZeros(K2_mat,8*3,8*3);
  femDoubleMat K3_mat;
  femUtils::matZeros(K3_mat,8*3,8*3);
  femDoubleMat K4_mat;
  femUtils::matZeros(K4_mat,8*3,8*3);

  // Shape functions, global derivatives and geometric matrix
  femDoubleVec shapeFunction;
  double detJ = 0.0;
  double tau_SUPG = 0.0;
  double nu_C = 0.0;
  double el_volume = 0.0;
  double gp_lam = 0.0;
  femDoubleMat shapeDeriv;
  femUtils::matZeros(shapeDeriv,8,3);
  femDoubleMat shapeDeriv_0;
  femUtils::matZeros(shapeDeriv_0,8,3);
  femDoubleMat elGeomMat;

  // Local element solution
  femDoubleVec el_sol_p(8*3,0.0);
  // Mass times previous solution
  femDoubleVec m_sol_p(8*3,0.0);
  femDoubleVec H_vec(8*3,0.0);
  femDoubleVec u_gp(3,0.0);
  
  // Lumped Bubnov Mass
  femDoubleVec Lumped_BMass(8*3,0.0);

  // Max and min velocity modules
  double max_vel_mod = 0.0;
  double min_vel_mod = 0.0;
  double min_vel_X = 0.0;
  double max_vel_X = 0.0;
  double min_vel_Y = 0.0;
  double max_vel_Y = 0.0;
  double min_vel_Z = 0.0;
  double max_vel_Z = 0.0;
  // Max and min tau coefficient
  double max_tau = 0.0;
  double min_tau = 0.0;
  // Max and min lamda value
  double min_lam = 0.0;
  double max_lam = 0.0;
  // Max and min E,nu 
  double cfl_E_min = 0.0;
  double cfl_E_max = 0.0;
  double cfl_nu_min = 0.0;
  double cfl_nu_max = 0.0;
  // Element characteristic lengths
  double el_len_X = 0.0;
  double el_len_Y = 0.0;
  double el_len_Z = 0.0;
  // Viscous CFL
  double v_CFL_X = 0.0;
  double v_CFL_Y = 0.0;
  double v_CFL_Z = 0.0;
  double v_CFL = 0.0;
  // Advective CFL
  double a_CFL_X = 0.0;
  double a_CFL_Y = 0.0;
  double a_CFL_Z = 0.0;
  double a_CFL = 0.0;
  // Maximum adv and visc CFL
  double v_CFL_max = 0.0;
  double a_CFL_max = 0.0;
  // Maximum structural CFL
  double el_cfl_E = 0.0;
  double el_cfl_nu = 0.0;
  // Max and min K matrix norms
  double min_K_norm = 0.0;
  double max_K_norm = 0.0;

  // ==============================================
  // PRELIMIANRY ASSEMBLY OF THE LUMPED BUBNOV MASS
  // ==============================================
  form_lumped_mass(model, rho, l_M);
  
  // INIT LINEAR SOLVER 
  femSolver* lin_solve;
  if(solve_static_K){  
    lin_solve = new femSolver();
  }
  
  // femUtils::writeVectorToFile("lumped_mass.txt",l_M);
  // getchar();

  // Get fluid properties
  // printf("%e %e %e %e\n",rho,mu,lam,cI);

  // ==============
  // RESET LOG FILE
  // ==============
  // LOG FILE POINTER
  FILE* log_fptr = fopen(model->log_file.c_str(), "w");
  fclose(log_fptr); 
  
  // =========
  // TIME LOOP
  // =========
  for(ulint loopTime=1;loopTime<model->totalSteps;loopTime++){    

    if((loopTime == 1) || ((loopTime % 10) == 0)){
      write_log_header(model->log_file);
    }

    // Update current time
    currentTime += model->timeStep;

    for(ulint loopA=0;loopA<sol_n.size();loopA++){
      sol_n[loopA][0] = 0.0;
      sol_n[loopA][1] = 0.0;
      sol_n[loopA][2] = 0.0;
    }    

    if(solve_static_K){  
      for(ulint loopA=0;loopA<totNodes*3;loopA++){
        glob_f[loopA] = 0.0;
        for(ulint loopB=0;loopB<totNodes*3;loopB++){
          glob_K[loopA][loopB] = 0.0;
        }
      }
    }

    // Initialize Lambda coefficient for each element
    for(ulint loopA=0;loopA<totElements;loopA++){
        if(use_constant_lamda){
            el_lam[loopA] = alpha;
        }else{
            el_lam[loopA] = 0.0;
        }
    }

    // Init max and min stabilization coefficient 
    // for current time steps    
    min_tau = std::numeric_limits<double>::max();
    max_tau = -std::numeric_limits<double>::max();

    // Init max and min of the manbda coefficient
    min_lam = std::numeric_limits<double>::max();
    max_lam = -std::numeric_limits<double>::max();

    // Init max and min CFL for current time step
    v_CFL_max = 0.0;
    a_CFL_max = 0.0;

    // Init min and max norms
    min_K_norm = std::numeric_limits<double>::max();
    max_K_norm = -std::numeric_limits<double>::max();

    // Init max and min E, nu
    cfl_E_min  =  std::numeric_limits<double>::max();
    cfl_E_max  = -std::numeric_limits<double>::max();
    cfl_nu_min =  std::numeric_limits<double>::max();
    cfl_nu_max = -std::numeric_limits<double>::max();
    
    // ============
    // ELEMENT LOOP
    // ============

    for(ulint loopElement=0;loopElement<model->elementList.size();loopElement++){

      // Get total number of nodes in element
      femIntVec el_conn = model->elementList[loopElement]->elementConnections;
      ulint tot_el_nodes = el_conn.size();

      // Get solution from previous time step
      for(ulint loopA=0;loopA<tot_el_nodes;loopA++){
        curr_node = model->elementList[loopElement]->elementConnections[loopA];
        el_sol_p[loopA*3 + 0] = sol_p[curr_node][0];
        el_sol_p[loopA*3 + 1] = sol_p[curr_node][1];
        el_sol_p[loopA*3 + 2] = sol_p[curr_node][2];
      }

      // printf("Element: %ld\n",loopElement);
      // femUtils::printVector(el_sol_p);
      // getchar();

      // SET ELEMENT-BASED SOLUTION AND MATRICES TO ZERO
      for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
        for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
          Mass_mat[loopA][loopB] = 0.0;
          K_mat[loopA][loopB] = 0.0;
          C_mat[loopA][loopB] = 0.0;
          // Contributions
          K1_mat[loopA][loopB] = 0.0;
          K2_mat[loopA][loopB] = 0.0;
          K3_mat[loopA][loopB] = 0.0;
          K4_mat[loopA][loopB] = 0.0;
        }
      }

      // GET CHARACTERISTIC LENGTHS IN X,Y,Z FOR CURRENT ELEMENT
      get_el_lengths(model,loopElement,el_len_X,el_len_Y,el_len_Z);

      //COMPUTE THE VISCOUS CFL
      v_CFL_X = 2.0*(mu/rho)*(model->timeStep/(el_len_X*el_len_X));
      v_CFL_Y = 2.0*(mu/rho)*(model->timeStep/(el_len_Y*el_len_Y));
      v_CFL_Z = 2.0*(mu/rho)*(model->timeStep/(el_len_Z*el_len_Z));
      v_CFL   = sqrt(v_CFL_X*v_CFL_X + v_CFL_Y*v_CFL_Y + v_CFL_Z*v_CFL_Z);
      if(v_CFL > v_CFL_max){
        v_CFL_max = v_CFL;
      }

      // COMPUTE AVERAGE B MATRIX
      if(use_B_hat){
        eval_avg_shapeDeriv(model,loopElement,shapeDeriv_0);
      }
      
      // Get Gauss Points and Weights
      femDoubleMat gps = rule->getCoords(tot_el_nodes,d3);
      femDoubleVec gpw = rule->getWeights(tot_el_nodes,d3);

      // Init element volume for integration
      el_volume = 0.0;

      // ================
      // GAUSS POINT LOOP
      // ================
      // printf("Total Gauss Points: %d\n",rule->getTotGP(tot_el_nodes,d3));  
      for(int igp=0;igp<rule->getTotGP(tot_el_nodes,d3);igp++){

        // printf("Gauss coord: %f %f %f\n",gps[igp][0],gps[igp][1],gps[igp][2]);
        // printf("Gauss weigth: %f\n",gpw[igp]);
        // getchar();

        // Eval Shape Function
        model->elementList[loopElement]->evalShapeFunction(model->nodeList,
                                                           gps[igp][0],gps[igp][1],gps[igp][2],
                                                           shapeFunction);

        // Eval velocity at current Gauss point
        u_gp[0] = 0.0; u_gp[1] = 0.0; u_gp[2] = 0.0;
        for(ulint loopA=0;loopA<tot_el_nodes;loopA++){
          u_gp[0] += shapeFunction[loopA]*el_sol_p[loopA*3 + 0];
          u_gp[1] += shapeFunction[loopA]*el_sol_p[loopA*3 + 1];
          u_gp[2] += shapeFunction[loopA]*el_sol_p[loopA*3 + 2];
        }

        // printf("Velocity at Gauss point %d, %f %f %f\n",igp,u_gp[0],u_gp[1],u_gp[2]);

        // Compute advective CFL at current Gauss point
        a_CFL_X = fabs(u_gp[0]*model->timeStep)/(el_len_X);
        a_CFL_Y = fabs(u_gp[1]*model->timeStep)/(el_len_Y);
        a_CFL_Z = fabs(u_gp[2]*model->timeStep)/(el_len_Z);
        a_CFL   = sqrt(a_CFL_X*a_CFL_X + a_CFL_Y*a_CFL_Y + a_CFL_Z*a_CFL_Z);
        if(a_CFL > a_CFL_max){
          a_CFL_max = a_CFL;
        }

        // Eval Current Shape Derivatives Matrix
        model->elementList[loopElement]->evalGlobalShapeFunctionDerivative(model->nodeList,
                                                                           gps[igp][0],gps[igp][1],gps[igp][2],
                                                                           detJ,shapeDeriv);

        //femUtils::printMatrix(shapeDeriv);
        //getchar();

        // Add current Gauss point contribution to the 
        el_volume += detJ * gpw[igp];


        // Eval Geometric Element Matrix
        model->elementList[loopElement]->evalGeometricMatrix(model->nodeList,
                                                             gps[igp][0],gps[igp][1],gps[igp][2],
                                                             elGeomMat);

        // femUtils::printMatrix(elGeomMat);
        // getchar();

        // COMPUTE STABILIZATION COEFFICIENTS
        eval_stabilization(u_gp,rho,mu,cI,model->timeStep,elGeomMat,tau_SUPG,nu_C);        
        if(tau_SUPG > max_tau){
          max_tau = tau_SUPG;
        }
        if(tau_SUPG < min_tau){
          min_tau = tau_SUPG;
        }
        // printf("stabilization: tau_SUPG: %f, nu_C: %f",tau_SUPG,nu_C);
        // getchar();

        // COMPUTE LOCAL BULK MODULUS
        if(!use_constant_lamda){
            eval_bulk_modulus(u_gp,rho,mu,alpha,model->timeStep,elGeomMat,
                              el_len_X,el_len_Y,el_len_Z,
                              gp_lam);
        }else{
            gp_lam = el_lam[loopElement];
        }
        if(gp_lam > max_lam){
          max_lam = gp_lam;
        }
        if(gp_lam < min_lam){
          min_lam = gp_lam;
        }
        // Add Integral of lam over the current element
        el_lam[loopElement] += gp_lam * detJ * gpw[igp];

        //printf("%e %e\n",tau_SUPG,nu_C);
        //getchar();

        // Eval L matrix
        eval_L_matrix(shapeFunction,L_mat);

        //femUtils::printMatrix(L_mat);
        //getchar();

        // Eval U matrix
        eval_U_matrix(u_gp,U_mat);

        //femUtils::printMatrix(U_mat);
        //getchar();
    
        // Eval S matrix
        eval_S_matrix(shapeDeriv,S_mat);

        //femUtils::printMatrix(S_mat);
        //getchar();

        // Conpute the US product: result is 3x24
        for(ulint loopA=0;loopA<3;loopA++){
          for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<9;loopC++){
              temp_1 += U_mat[loopA][loopC] * S_mat[loopC][loopB];
            }
            US_mat[loopA][loopB] = temp_1;
          }
        }

        //femUtils::printMatrix(US_mat);
        //getchar();

        // Eval B matrix
        eval_B_matrix(use_B_hat,shapeDeriv,shapeDeriv_0,B_mat);

        //femUtils::printMatrix(B_mat);
        //getchar();

        // Eval D matrix
        eval_D_matrix(gp_lam,mu,D_mat);

        //femUtils::printMatrix(D_mat);
        //getchar();

        // Conpute DB Matrix
        for(ulint loopA=0;loopA<6;loopA++){
          for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<6;loopC++){
              temp_1 += D_mat[loopA][loopC] * B_mat[loopC][loopB];
            }
            DB_mat[loopA][loopB] = temp_1;
          }
        }

        // Eval H matrix
        eval_H_matrix(shapeDeriv,H_vec);

        // FOR MASS, STIFFNESS AND STABILIZATION MATRICES
        for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
          for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
            
            // ===========
            // MASS MATRIX
            // ===========
            
            // BUBNOV MASS MATRIX
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<3;loopC++){
              temp_1 += rho * L_mat[loopC][loopA] * L_mat[loopC][loopB] * detJ * gpw[igp];
            }
            
            // STABILIZATION MASS MATRIX
            temp_2 = 0.0;
            for(ulint loopC=0;loopC<3;loopC++){
              temp_2 += tau_SUPG * US_mat[loopC][loopA] * L_mat[loopC][loopB] * detJ * gpw[igp];
            }

            // NO STABILIZATION AND NO CONSISTENT MASS ON RHS
            // TO AVOID DIFFUSION-LIKE COUPLING
            Mass_mat[loopA][loopB] += temp_1;

            // ================
            // STIFFNESS MATRIX
            // ================

            // BUBNOV STIFFNESS
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<3;loopC++){
              temp_1 += rho * L_mat[loopC][loopA] * US_mat[loopC][loopB] * detJ * gpw[igp];
            }
            K1_mat[loopA][loopB] += temp_1;

            // INCOMPRESSIBLE ELASTIC-TYPE STIFFNESS
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<6;loopC++){
              temp_1 += B_mat[loopC][loopA] * DB_mat[loopC][loopB] * detJ * gpw[igp];
            }
            K2_mat[loopA][loopB] += temp_1;

            // DOUBLE CONVECTIVE STIFFNESS
            temp_1 = 0.0;
            for(ulint loopC=0;loopC<3;loopC++){
              temp_1 += rho * tau_SUPG * US_mat[loopC][loopA] * US_mat[loopC][loopB] * detJ * gpw[igp];
            }
            K3_mat[loopA][loopB] += temp_1;

            // Final term with nu_C
            K4_mat[loopA][loopB] += rho * nu_C * H_vec[loopA] * H_vec[loopB] * detJ * gpw[igp];

          }
        } // END OF MATRIX COMPONENT LOOP

        //femUtils::printMatrix(Mass_mat);
        //getchar();

        //femUtils::printMatrix(K_mat);
        //getchar();

      } // END OF GAUSS POINT LOOP

      // getchar();

      // ADD COMPONENTS OF STIFFNESS MATRIX
      for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
          for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
              K_mat[loopA][loopB] = 0.0;
              if(include_K1){
                // BUBNOV CONVECTION
                K_mat[loopA][loopB] += K1_mat[loopA][loopB];
              }
              if(include_K2){
                // DIFFUSION BDB
                K_mat[loopA][loopB] += K2_mat[loopA][loopB];
              }
              if(include_K3){
                // DOUBLE CONVECTIVE STABILIZATION
                K_mat[loopA][loopB] += K3_mat[loopA][loopB];
              }
              if(include_K4){
                // DIVERGENCE STABILIZATION
                K_mat[loopA][loopB] += K4_mat[loopA][loopB];
              }
          }
      }

      //printf("MASS\n");
      //femUtils::printMatrix(Mass_mat);
      //getchar();

      //printf("K\n");
      //femUtils::printMatrix(K_mat);
      //getchar();

      if(femUtils::getMatrixNorm(K_mat) > max_K_norm){
          max_K_norm = femUtils::getMatrixNorm(K_mat);
      }
      if(femUtils::getMatrixNorm(K_mat) < min_K_norm){
          min_K_norm = femUtils::getMatrixNorm(K_mat);
      }

      // printf("Element: %ld\n",loopElement);
      // printf("K1 norm: %e\n",);
      // printf("K2 norm: %e\n",femUtils::getMatrixNorm(K2_mat));
      // getchar();
      // printf("K3 norm: %e\n",femUtils::getMatrixNorm(K3_mat));
      // printf("K4 norm: %e\n",femUtils::getMatrixNorm(K4_mat));
      // getchar();

      // printf("Mass matrix norm: %f\n",femUtils::getMatrixNorm(Mass_mat));
      // printf("K matrix norm: %f\n",femUtils::getMatrixNorm(K_mat));
      // getchar();

      // GET LUMPED MASS FROM FULLY ASSEMBLE MATRIX
      // for(ulint loopA=0;loopA<tot_el_nodes;loopA++){
      //   curr_node = model->elementList[loopElement]->elementConnections[loopA];
      //   for(ulint loopB=0;loopB<3;loopB++){
      //     Lumped_BMass[loopA*3 + loopB] = l_M[curr_node*3 + loopB];
      //   }
      // }

      // Normalize the integral of lam over the element
      el_lam[loopElement] = (el_lam[loopElement]/el_volume);

      // COMPUTE THE STRUCTURAL EQUIVALENT E,nu FOR THE CURRENT ELEMENT
      el_cfl_E = (mu*(3*el_lam[loopElement] + 2*mu))/(el_lam[loopElement] + mu);
      el_cfl_nu = el_lam[loopElement]/(2*(el_lam[loopElement] + mu));
      if(el_cfl_E < cfl_E_min){
          cfl_E_min = el_cfl_E;
      }
      if(el_cfl_E > cfl_E_max){
          cfl_E_max = el_cfl_E;
      }
      if(el_cfl_nu < cfl_nu_min){
          cfl_nu_min = el_cfl_nu;
      }
      if(el_cfl_nu > cfl_nu_max){
          cfl_nu_max = el_cfl_nu;
      }

      // ASSEMBLE THE FORCE VECTORS
      for(ulint loopA=0;loopA<tot_el_nodes*3;loopA++){
        k_sol_p[loopA][loopElement] = 0.0;
        for(ulint loopB=0;loopB<tot_el_nodes*3;loopB++){
          // Careful to the sign
          k_sol_p[loopA][loopElement] += - model->timeStep * K_mat[loopA][loopB] * el_sol_p[loopB];
        }
      }

      // ADD TO GLOBAL SOLUTION
      for(ulint loopNode=0;loopNode<model->elementList[loopElement]->elementConnections.size();loopNode++){

          curr_node = model->elementList[loopElement]->elementConnections[loopNode];

          for(ulint loopDOF=0;loopDOF<3;loopDOF++){

              if(time_alg == eaForwardEuler){

                sol_n[curr_node][loopDOF] += k_sol_p[loopNode*3 + loopDOF][loopElement];

              }else if(time_alg == eaAdamsBashforth){

                if(loopTime == 1){

                  sol_n[curr_node][loopDOF] += k_sol_p[loopNode*3 + loopDOF][loopElement];

                }else{

                  sol_n[curr_node][loopDOF] += (3.0/2.0)*k_sol_p[loopNode*3 + loopDOF][loopElement] - (1.0/2.0)*k_sol_pp[loopNode*3 + loopDOF][loopElement];

                }
              }              
          }    
      }

      // getchar();

      // TEST - ASSEMBLE IN GLOBAL STIFFNESS
      if(solve_static_K){
        assemble_in_glob(model,loopElement,K_mat,glob_K);
      }  
      
    } // END OF ELEMENT LOOP    

    // // UPDATE SOLUTION VECTOR
    for(ulint loopA=0;loopA<sol_n.size();loopA++){
      for(ulint loopB=0;loopB<3;loopB++){
        sol_n[loopA][loopB] = sol_p[loopA][loopB] + (1.0/l_M[loopA*3 + loopB]) * sol_n[loopA][loopB];
      }
    }

    // SET WALL BOUNDARY CONDITIONS
    model->setDirichletBC(sol_n);
    // SET BOUNDARY CONDITIONS FOR SOL_N
    model->setNodeVelocity(currentTime,sol_n);

    // POST-PROCESS PRESSURE
    compute_pressure(model,el_lam,sol_n,pres_n);

    // SET GLOBALLY FIXED DOFs
    for(ulint loopA=0;loopA<sol_p.size();loopA++){
        if(global_fixed_dof != "NONE"){
            if(global_fixed_dof == "X"){

                sol_n[loopA][0] = 0.0;

            }else if(global_fixed_dof == "Y"){

                sol_n[loopA][1] = 0.0;

            }else if(global_fixed_dof == "Z"){

                sol_n[loopA][2] = 0.0;

            }else if(global_fixed_dof == "XY"){
            
                sol_n[loopA][0] = 0.0;
                sol_n[loopA][1] = 0.0;

            }else if(global_fixed_dof == "YZ"){

                sol_n[loopA][1] = 0.0;
                sol_n[loopA][2] = 0.0;

            }else if(global_fixed_dof == "XZ"){

                sol_n[loopA][0] = 0.0;
                sol_n[loopA][2] = 0.0;

            }else{
              printf("ERROR: Invalid global fixed DOFs string.\n");
              exit(-1);
            }
        }
    }

    // APPLY DIRICHLET CONDITIONS TO MATRIX
    if(solve_static_K){
      printf("Applying Dirichlet Conditions...");
      fflush(stdout);
      apply_dir_vel_glob_K_f(model,currentTime,glob_K,glob_f);
      printf("OK\n");
      fflush(stdout);

      // SOLVE SYSTEM
      printf("Solving linear system...");
      fflush(stdout);
      lin_sol = lin_solve->solveLinearSystem(totNodes*3,totNodes*3, glob_K,glob_f);
      printf("OK\n");
      fflush(stdout);
    }

    // UPDATE PREV SOLUTION
    for(ulint loopA=0;loopA<sol_p.size();loopA++){
      sol_p[loopA][0] = sol_n[loopA][0];
      sol_p[loopA][1] = sol_n[loopA][1];
      sol_p[loopA][2] = sol_n[loopA][2];
    }

    // UPDATE FORCE VECTORS FOR PREVIOUS TIME STEPS
    for(ulint loopA=0;loopA<totElements;loopA++){
      for(ulint loopB=0;loopB<8*3;loopB++){
        k_sol_pp[loopB][loopA] = k_sol_p[loopB][loopA];
      }
    }

    // Compute current solution extrema
    get_solution_extrema(sol_n,
                         min_vel_mod,max_vel_mod,
                         min_vel_X,max_vel_X,
                         min_vel_Y,max_vel_Y,
                         min_vel_Z,max_vel_Z);

    // PRINT ITERATION MESSAGE
    printf("%5ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3e %10.3e %10.3e %10.3e\n",loopTime,
                                                                                                  currentTime,
                                                                                                  min_vel_X,max_vel_X,
                                                                                                  min_vel_Y,max_vel_Y,
                                                                                                  min_vel_Z,max_vel_Z,
                                                                                                  min_tau,max_tau,
                                                                                                  min_lam,max_lam,
                                                                                                  v_CFL_max,
                                                                                                  a_CFL_max,
                                                                                                  cfl_E_min,cfl_E_max,
                                                                                                  cfl_nu_min,cfl_nu_max);
    fflush(stdout);
    FILE* fptr = fopen(model->log_file.c_str(), "a");
    
    // Write some text to the file
    fprintf(fptr, "%5ld %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3f %10.3f %10.3e %10.3e %10.3e %10.3e\n",loopTime,
                                                                                                  currentTime,
                                                                                                  min_vel_X,max_vel_X,
                                                                                                  min_vel_Y,max_vel_Y,
                                                                                                  min_vel_Z,max_vel_Z,
                                                                                                  min_tau,max_tau,
                                                                                                  min_lam,max_lam,
                                                                                                  v_CFL_max,
                                                                                                  a_CFL_max,
                                                                                                  cfl_E_min,cfl_E_max,
                                                                                                  cfl_nu_min,cfl_nu_max);
    // Close the file
    fclose(fptr); 

    // getchar();

    // EXIT CONDITION ON THE MAX VELOCITY
    if(max_vel_mod > max_allowed_vel){
      FILE* fptr = fopen(model->log_file.c_str(), "a");
      fprintf(fptr,"ERROR: Maximum velocity module larger than admissible threshold.\n");
      fprintf(fptr,"SOLUTION TERMINATED.\n");
      fclose(fptr); 
      printf("ERROR: Maximum velocity module larger than admissible threshold.\n");
      printf("SOLUTION TERMINATED.\n");
      exit(-1);
    }
  
    // printf("Max and min K1 norms: min: %e, max: %e\n",min_K_norm,max_K_norm);
           
    // SAVE RESULTS IF REQUESTED
    if(saveCounter == model->saveEvery){      

      printf("Saving results at time step: %ld\n",loopTime);
      fflush(stdout);

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
        res_temp.clear();
        res_temp.push_back(sol_n[loopA][0]);
        res_temp.push_back(sol_n[loopA][1]);
        res_temp.push_back(sol_n[loopA][2]);
        res->values.push_back(res_temp);
      }
      model->resultList.push_back(res);

      // SAVE PRESSURE
      res = new femResult();
      res->label = string("pressure");
      res->type = frElement;
      res->numComponents = 1;
      // Assign to values
      for(size_t loopA=0;loopA<totElements;loopA++){
        res_temp.clear();
        res_temp.push_back(pres_n[loopA]);
        res->values.push_back(res_temp);
      }
      model->resultList.push_back(res);

      // SAVE STATIC SOLUTION
      if(solve_static_K){
        res = new femResult();
        res->label = string("static_sol");
        res->type = frNode;
        res->numComponents = 3;
        // Assign to values
        for(size_t loopA=0;loopA<totNodes;loopA++){
          res_temp.clear();
          res_temp.push_back(lin_sol[loopA*3 + 0]);
          res_temp.push_back(lin_sol[loopA*3 + 1]);
          res_temp.push_back(lin_sol[loopA*3 + 2]);
          res->values.push_back(res_temp);
        }
        model->resultList.push_back(res);
      }
      // EXPORT MODEL
      model->ExportToVTKLegacy(string(model->sol_file_prefix + "_" + femUtils::intToStr(loopTime) + ".vtk"),false);
    }
    // Update counter
    saveCounter++;
  } // END OF TIME STEP LOOP
}




