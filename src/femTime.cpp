#include "femTime.h"

double femTime::LocalSFTime = 0.0;
double femTime::JacobianMatTime = 0.0;
double femTime::JacInversionTime = 0.0;
double femTime::globShDerivAllocTime = 0.0;
double femTime::totLHSAssemblyTime = 0.0;
double femTime::totLHSMatAssembleTime = 0.0;
double femTime::totTrilinosCompleteFillTime = 0.0;
double femTime::multShDerivAllocTime = 0.0;
double femTime::writeLHSAssemblyProgressTime = 0.0;
double femTime::totAssembleFromElementList = 0.0;
double femTime::evalGlobalShapeFunctionDerivativeTime = 0.0;
double femTime::evalLocalMatrixTime = 0.0;

// Print the Timem infos to screen
void femTime::printToScreen(){
  printf("Local Shape Function Construction Time: %f [s]\n",LocalSFTime);
  printf("Jacobian Matrix Construction Time: %f [s]\n",JacobianMatTime);
  printf("Jacobian Matrix Inversion Time: %f [s]\n",JacInversionTime);
  printf("globShDeriv Allocation Time: %f [s]\n",globShDerivAllocTime);
  printf("multShDeriv Allocation Time: %f [s]\n",multShDerivAllocTime);
  printf("LHS Matrix Assemble Time: %f [s]\n",totLHSMatAssembleTime);
  printf("Trilinos Complete Fill Time: %f [s]\n",totTrilinosCompleteFillTime);
  printf("Write LHS Assembly Progress Time: %f [s]\n",writeLHSAssemblyProgressTime);
  printf("Total Assembly From Element List Time: %f [s]\n",totAssembleFromElementList);
  printf("Total evalGlobalShapeFunctionDerivativeTime Time: %f [s]\n",evalGlobalShapeFunctionDerivativeTime);
  printf("Total evalLocalMatrixTime Time: %f [s]\n",evalLocalMatrixTime);
  printf("--- Total LHS Assembly Time: %f [s]\n",totLHSAssemblyTime);
}
