#include "femCMMExplicitSolidSolver.h"

const double w[] = {1.0/3.0, 1.0/3.0, 1.0/3.0};
const double lN[3][3] = {{0.5, 0.0, 0.5},
                         {0.0, 0.5, 0.5},
                         {0.5, 0.5, 0.0}};
const double lDN[3][4] = {{-1.0, 1.0, 0.0, 0.0},
                          {-1.0, 0.0, 1.0, 0.0},
                          {-1.0, 0.0, 0.0, 1.0}};

// CONSTRUCTOR
femCMMExplicitSolidSolver::femCMMExplicitSolidSolver(){
}


femDoubleMat get_transformation(ulint iElm,
                                femModel* model){

  // COMPLETE!!!
  // femDoubleMat coords;
  // femDoubleMat mat;
  // femUtils::matZeros(mat,3,3);

  femDoubleMat transf;

  // // Read coords from 
  // coords = get_elm_coords(model);

  // // Get first axis
  // axis_1[0] = coords[1][0] - coords[0][0];
  // axis_1[1] = coords[1][1] - coords[0][1];
  // axis_1[2] = coords[1][2] - coords[0][2];

  // // Get aux axis 
  // aux[0] = coords[2] - coords[1];
  // aux[1] = coords[2] - coords[1];
  // aux[2] = coords[2] - coords[1];

  // Compute normal

  return transf;
}

femDoubleMat get_constitutive_tensor(double E,double nu, double k){
  femDoubleMat res;
  femUtils::matZeros(res,5,5);
  // Fill tensor
  res[0][0] = 1.0;
  res[0][1] = nu;
  //
  res[1][0] = nu;
  res[1][1] = 1.0;
  //
  res[2][2] = 0.5*(1.0-nu);
  //
  res[3][3] = 0.5*k*(1.0-nu);
  //
  res[4][4] = 0.5*k*(1.0-nu);
  // Return tensor
  return res;
}

// Compute residual for current element and assemble in node-wise array
void assemble_RES(ulint iElm,
                  femModel* model,
                  // Displacement solution
                  const femDoubleMat& disps,
                  const femDoubleVec& el_vol,
                  // Return nodal residual
                  femDoubleMat& RES){

  // double wGp = 0.0;
  // double temp = 0.0;
  // double eps[3][3];
  // double sigma[3][3];
  // double lRes[3][3];
  // ulint curr_node = 0;
  // // Assemble structural residual for current element

  // // Get constitutive tensor
  // double young   = model->CMMProps[0];
  // double poisson = model->CMMProps[1];
  // double k       = model->CMMProps[2];

  // // Get constitutive and transformation matrix
  // femDoubleMat constitutive = get_constitutive_tensor(young,poisson,k);
  // femDoubleMat trans_mat    = get_transformation(iElm,model);

  // // Initialize local residual
  // for(ulint i=0;i<3;i++){
  //   for(ulint j=0;j<3;j++){
  //     lRes[i][j] = 0.0;
  //   }
  // }

  // // Get element volume
  // double Ve = el_vol[iElm];

  // // Gauss point loop
  // for(ulint iGp=0;iGp<4;iGp++){

  //   // Get Gauss point weight
  //   wGp = w[iGp] * Ve;

  //   // Evaluate stress tensor at current Gauss point
  //   // Evaluate infinitesimal strains at iGp
  //   for(ulint i=0;i<3;i++){
  //     for(ulint j=0;j<3;j++){
  //       eps[i][j] = 
  //     }
  //   }

  //   // Compute local stress tensor
  //   for(ulint loopA=0;loopA<3;loopA++){
  //     for(ulint loopB=0;loopB<3;loopB++){
  //       temp = 0.0;
  //       for(ulint loopC=0;loopC<3;loopC++){
  //         temp += constitutive[loopA][loopC] * eps[loopC][loopB];
  //       }
  //       sigma[loopA][loopB] = temp;
  //     }
  //   }

  //   // Rotate tensor to global Cartesian system

  //   trans_mat[][] * sigma[loopA][loopB] * trans_mat[][];

  //   // Compute divergence

  //   for(ulint iNode=0;iNode<3;iNode++){
  //     for(ulint j=0;j<3;j++){
  //       // Assemble internal forces
  //       lRes[iNode][j] += wGp*DivSigma[iNode][j];
  //       // Assemble external distributed forces
  //       lRes[iNode][j] += wGp*force[iNode][j]*lN[][];
  //     }
  //   }
  // }

  // // Assemble the local residual to the nodes
  // for(ulint i=0;i<3;i++){
  //   curr_node = model->elementList[iElm]->elementConnections[i];
  //   for(ulint j=0;j<3;j++){
  //     RES[curr_node][j] += lRes[i][j];
  //   }
  // }
}

// EXPLICIT CMM SOLVER
void femCMMExplicitSolidSolver::solve(femModel* model){

  // double currentTime = 0.0;
  // // Totals
  // double totElements = model->elementList.size();
  // double totNodes = model->nodeList.size();

  // // Declare and initialize variables
  // femDoubleMat sol_v;
  // femUtils::matZeros(sol_v,totNodes,3);
  // femDoubleMat sol_v_prev;
  // femUtils::matZeros(sol_v_prev,totNodes,3);
  // femDoubleMat sol_disp;
  // femUtils::matZeros(sol_disp,totNodes,3);
  // femDoubleMat residual;
  // femUtils::matZeros(residual,totNodes,3);
  // femDoubleVec mass(totNodes,0.0);

  // // Apply initial velocity/displacement

  // // TIME LOOP
  // for(uint loopTime=1;loopTime<model->totalSteps;loopTime++){

  //   // Update current time
  //   currentTime += model->timeStep;

  //   // Set residual to zero
  //   for(ulint loopA=0;loopA<totNodes;loopA++){
  //     for(ulint loopB=0;loopB<3;loopB++){
  //       residual[loopA][loopB] = 0.0;
  //     }
  //   }
    

  //   // Element Loop
  //   for(ulint loopElement=0;loopElement<totElements;loopElement++){


  //     // Assemble residual for the current element
  //     assemble_RES(residual);

  //   }

  //   // Update quantities
  //   for(ulint loopA=0;loopA<totNodes;loopA++){
  //     for(ulint j=0;j<3;j++){
  //       // Update nodal velocity
  //        sol_v[loopA][j] = sol_v_prev[loopA][j] - model->timeStep*residual[loopA][j]/mass[loopA];
  //        // Compute diplacements 
  //        sol_disp[loopA][j] = (sol_v[loopA][j] - sol_v_prev[loopA][j])/model->timeStep;
  //     }
  //   }

  //   // Apply boundary conditions



    

  // }
}



