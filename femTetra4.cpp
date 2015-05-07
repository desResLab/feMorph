# include "femConstants.h"
# include "femUtils.h"
# include "femElement.h"

// ====================================
// Eval Tetra4 Volume from Coord Matrix
// ====================================
double EvalTetVolumeFromCoordMat(double** coordMat){
  // Allocate
  double firstVec[3] = {0.0};
  double secondVec[3] = {0.0};
  double thirdVec[3] = {0.0};

  // Get The Vectors
  // First
  firstVec[0] = coordMat[0][1] - coordMat[0][0];
  firstVec[1] = coordMat[1][1] - coordMat[1][0];
  firstVec[2] = coordMat[2][1] - coordMat[2][0];
  // Second
  secondVec[0] = coordMat[0][2] - coordMat[0][0];
  secondVec[1] = coordMat[1][2] - coordMat[1][0];
  secondVec[2] = coordMat[2][2] - coordMat[2][0];
  // Third
  thirdVec[0] = coordMat[0][3] - coordMat[0][0];
  thirdVec[1] = coordMat[1][3] - coordMat[1][0];
  thirdVec[2] = coordMat[2][3] - coordMat[2][0];

  // Eval External Product
  double auxVec[3] = {0.0};
  femUtils::Do3DExternalProduct(secondVec,thirdVec,auxVec);

  // Eval Internal Product
  return (femUtils::Do3DInternalProduct(firstVec,auxVec)/6.0);
}

// ========================
// Get External Tet Volumes
// ========================
void EvalExternalTetVolumes(double* pointCoords, double** coordMat, double* extTetVol){
  // Allocate current coordinate matrix
  double** currCoordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    currCoordMat[loopA] = new double[kTetra4Nodes];
  }
  // Eval Volume Coords
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    for(int loopB=0;loopB<kTetra4Nodes;loopB++){
      if (loopB != loopA){
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = coordMat[loopC][loopB];
        }
      }else{
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = pointCoords[loopC];
        }
      }
    }
    // Eval Volume Using CurrCoordMat
    extTetVol[loopA] = EvalTetVolumeFromCoordMat(currCoordMat);
  }
  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] currCoordMat[loopA];
  }
  delete [] currCoordMat;
}

// =======================================
// Assemble Tetra4 Coordinates in a Matrix
// =======================================
void femTetra4::AssembleTetCoordsMat(double dispFactor, std::vector<femNode*> &nodeList, double** coordMat){
  // Fill Matrix
  int currNode = 0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    // Get Current Node
    currNode = elementConnections[loopA];
    // Store In Matrix
    for(int loopB=0;loopB<3;loopB++){
      coordMat[loopB][loopA] = nodeList[currNode]->coords[loopB] + dispFactor * nodeList[currNode]->displacements[loopB];
    }
  }
}

// ================================
// ASSEMBLE STABILIZATION PARAMETER
// ================================
/*virtual void femTetra4::assembleStab(femOptions options, femIntegrationRule* rule,femDoubleMat &nodeVelocities, std::vector<double> &tauSUPG){
  // Loop Through the Gauss Points
  for(int loopGauss=0;loopGauss<rule->totalGP;loopGauss++){
    evalShapeFunctionDerivative(nodeList,rule->coords[0],rule->coords[1],rule->coords[2],shapeDeriv);
    // Assemble First Term
    firstTerm[loopGauss] = 0.0;
    for(int loopA=0;loopA<numberOfNodes;loopA++){
      dotProduct = 0.0;
      currNode = elementConnections[loopA];
      for(int loopB=0;loopB<kDims;loopB++){
        dotProduct += nodeVelocities[currNode][loopB] * shapeDeriv[loopA][loopB];
      }
      firstTerm[loopGauss] += fabs(dotProduct)
    }
    firstTerm[loopGauss] = (1.0/(firstTerm[loopGauss]));
    // Assemble Second Term
    secondTerm[loopGauss] = (options.deltaT/2.0);

  }
}*/

// ===================
// Eval Element Volume
// ===================
double femTetra4::EvalVolume(double dispFactor, std::vector<femNode*> &nodeList){
  // Put The Coordinates in a Matrix
    // Allocate current coordinate matrix
  double** coordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    coordMat[loopA] = new double[kTetra4Nodes];
  }
  AssembleTetCoordsMat(dispFactor,nodeList,coordMat);

  // Eval Tethrahedral Volume
  return EvalTetVolumeFromCoordMat(coordMat);
}

// =================================================
// SWAP ELEMENT NODES TO RESTORE A POSITIVE JACOBIAN
// =================================================
void femTetra4::fixConnectivities(std::vector<femNode*> &nodeList){
  int temp = elementConnections[1];
  elementConnections[1] = elementConnections[2];
  elementConnections[2] = temp;
}

// ========================================
// Eval Volume Coordinates of Point in Tet4
// ========================================
void femTetra4::EvalVolumeCoordinates(double dispFactor, double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){

  // Put The Coordinates in a Matrix
    // Allocate current coordinate matrix
  double** coordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    coordMat[loopA] = new double[kTetra4Nodes];
  }
  AssembleTetCoordsMat(dispFactor,nodeList,coordMat);

  // Eval Tethrahedral Volume
  double TetVolume = EvalTetVolumeFromCoordMat(coordMat);

  // Eval Other Volumes
  double currVolumes[kTetra4Nodes] = {0.0};
  EvalExternalTetVolumes(pointCoords,coordMat,currVolumes);

  // Eval Shape Function Summation
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++)
  {
    sum += currVolumes[loopA];
  }
  //if (fabs(sum-1.0)>kMathZero){
  //  throw femException("Internal: Tet4 Volume Coordinates do not sum up to one.");
  //}

  // Compute Final Volume Coordinates
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    volCoords[loopA] = (currVolumes[loopA]/TetVolume);
  }

  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] coordMat[loopA];
  }
  delete [] coordMat;
}


// Check if Node is Inside
bool femTetra4::isNodeInsideElement(double dispFactor, double* pointCoords,std::vector<femNode*> &nodeList){
  // Init Result
  bool isInside = true;

  // Eval Volume Coordinates
  double volCoords[kTetra4Nodes] = {0.0};
  EvalVolumeCoordinates(dispFactor,pointCoords,nodeList,volCoords);

  // Compute Result
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    isInside = (isInside)&&(volCoords[loopA] >= 0.0);
  }

  // Check if the volume coordinates sum up to one
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    sum += volCoords[loopA];
  }
  // NO: I use it to evaluate if a node is inside
  //if (fabs(sum-1.0)>kMathZero){
  //  throw femException("Error: Internal. Tet4 Shape Function don't sum up to one");
  //}
  // Return
  return isInside;
}

// ==========================================
// FINITE ELEMENT ROUTINES FOR TETRA4 ELEMENT
// ==========================================
void femTetra4::evalShapeFunction(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3, femDoubleVec &shapeFunction){
  shapeFunction.clear();
  shapeFunction.push_back(coord1);
  shapeFunction.push_back(coord2);
  shapeFunction.push_back(coord3);
  shapeFunction.push_back(1.0-coord1-coord2-coord3);
}

// DERIVATIVES OF SHAPE FUNCTIONS RESPECT TO LOCAL COORDINATES
void femTetra4::evalLocalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,
                                                 femDoubleMat &shapeDeriv){
  //        nodi - dim
  femDoubleVec temp;
  temp.push_back(1.0);
  temp.push_back(0.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
  temp.push_back(0.0);
  temp.push_back(1.0);
  temp.push_back(0.0);
  shapeDeriv.push_back(temp);
  temp.push_back(0.0);
  temp.push_back(0.0);
  temp.push_back(1.0);
  shapeDeriv.push_back(temp);
  temp.push_back(-1.0);
  temp.push_back(-1.0);
  temp.push_back(-1.0);
  shapeDeriv.push_back(temp);
}


// NOT GENERAL
/*
// EVAL SHAPE FUNCTION DERIVATIVES
void femTetra4::evalShapeFunctionDerivative(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3,
                                            femDoubleMat &shapeDeriv){

  // Eval Coordinate Differences
  double x12 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x13 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[2]]->coords[0];
  double x14 = nodeList[elementConnections[0]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x21 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x24 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x31 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x32 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x34 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[3]]->coords[0];
  double x42 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x43 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[2]]->coords[0];

  double y12 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y13 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[2]]->coords[1];
  double y14 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y21 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y23 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[2]]->coords[1];
  double y24 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y31 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[0]]->coords[1];
  double y32 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y34 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[3]]->coords[1];
  double y42 = nodeList[elementConnections[3]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y43 = nodeList[elementConnections[3]]->coords[1] - nodeList[elementConnections[2]]->coords[1];

  double z12 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z13 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[2]]->coords[2];
  double z14 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z21 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  double z23 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[2]]->coords[2];
  double z24 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z31 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[0]]->coords[2];
  double z32 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z34 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[3]]->coords[2];
  double z42 = nodeList[elementConnections[3]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z43 = nodeList[elementConnections[3]]->coords[2] - nodeList[elementConnections[2]]->coords[2];

  // Eval Jacobian Determinant
  double jac = (x21 * (y23 * z34 - y34 * z23) + x32 * (y34 * z12 - y12 * z34) + x43 * (y12 * z23 - y23 * z12));

  // Eval Global Derivatives of Volume Coordinates
  double globShapeDeriv[4][3];
  globShapeDeriv[0][0] = (y42*z32-y32*z42)/jac;
  globShapeDeriv[1][0] = (y31*z43-y34*z13)/jac;
  globShapeDeriv[2][0] = (y24*z14-y14*z24)/jac;
  globShapeDeriv[3][0] = (y13*z21-y12*z31)/jac;

  globShapeDeriv[0][1] = (x32*z42-x42*z32)/jac;
  globShapeDeriv[1][1] = (x43*z31-x13*z34)/jac;
  globShapeDeriv[2][1] = (x14*z24-x24*z14)/jac;
  globShapeDeriv[3][1] = (x21*z13-x31*z12)/jac;

  globShapeDeriv[0][2] = (x42*y32-x32*y42)/jac;
  globShapeDeriv[1][2] = (x31*y43-x34*y13)/jac;
  globShapeDeriv[2][2] = (x24*y14-x14*y24)/jac;
  globShapeDeriv[3][2] = (x13*y21-x12*y31)/jac;

  // Eval Shape function Derivatives respect to volume coordinate
  double locShapeDeriv[4][4];
  locShapeDeriv[0][0] = 1.0;
  locShapeDeriv[1][0] = 0.0;
  locShapeDeriv[2][0] = 0.0;
  locShapeDeriv[3][0] = -1.0;

  locShapeDeriv[0][1] = 0.0;
  locShapeDeriv[1][1] = 1.0;
  locShapeDeriv[2][1] = 0.0;
  locShapeDeriv[3][1] = -1.0;

  locShapeDeriv[0][2] = 0.0;
  locShapeDeriv[1][2] = 0.0;
  locShapeDeriv[2][2] = 1.0;
  locShapeDeriv[3][2] = -1.0;

  locShapeDeriv[0][3] = 0.0;
  locShapeDeriv[1][3] = 0.0;
  locShapeDeriv[2][3] = 0.0;
  locShapeDeriv[3][3] = 0.0;

  // Get Final Result
  // Loop on the Shape Function whose derivative needs to be evaluated
  for(int loopA=0;loopA<numberOfNodes;loopA++){
    // Loop on the global direction
    for(int loopB=0;loopB<kDims;loopB++){
      shapeDeriv[loopA][loopB] = 0.0;
      for(int loopC=0;loopC<numberOfNodes;loopC++){
        shapeDeriv[loopA][loopB] += globShapeDeriv[loopC][loopB]*locShapeDeriv[loopA][loopC];
      }
    }
  }

}
*/

// NOT GENERAL
/*double femTetra4::evalJacobian(std::vector<femNode*> nodeList, double coord1, double coord2, double coord3){
  // Define Quantities
  double x21 = nodeList[elementConnections[1]]->coords[0] - nodeList[elementConnections[0]]->coords[0];
  double x32 = nodeList[elementConnections[2]]->coords[0] - nodeList[elementConnections[1]]->coords[0];
  double x43 = nodeList[elementConnections[3]]->coords[0] - nodeList[elementConnections[2]]->coords[0];

  double y12 = nodeList[elementConnections[0]]->coords[1] - nodeList[elementConnections[1]]->coords[1];
  double y23 = nodeList[elementConnections[1]]->coords[1] - nodeList[elementConnections[2]]->coords[1];
  double y34 = nodeList[elementConnections[2]]->coords[1] - nodeList[elementConnections[3]]->coords[1];

  double z12 = nodeList[elementConnections[0]]->coords[2] - nodeList[elementConnections[1]]->coords[2];
  double z23 = nodeList[elementConnections[1]]->coords[2] - nodeList[elementConnections[2]]->coords[2];
  double z34 = nodeList[elementConnections[2]]->coords[2] - nodeList[elementConnections[3]]->coords[2];

  return (x21 * (y23 * z34 - y34 * z23) + x32 * (y34 * z12 - y12 * z34) + x43 * (y12 * z23 - y23 * z12));

}*/


