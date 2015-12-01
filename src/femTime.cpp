#include "femTime.h"

double femTime::LocalSFTime = 0.0;
double femTime::JacobianMatTime = 0.0;
double femTime::JacInversionTime = 0.0;
double femTime::globShDerivAllocTime = 0.0;
double femTime::totLHSAssemblyTime = 0.0;


// Print the Timem infos to screen
void femTime::printToScreen(){
  printf("Local Shape Function Construction Time: %f [s]\n",LocalSFTime);
  printf("Jacobian Matrix Construction Time: %f [s]\n",JacobianMatTime);
  printf("Jacobian Matrix Inversion Time: %f [s]\n",JacInversionTime);
  printf("globShDeriv Allocation Time: %f [s]\n",globShDerivAllocTime);
  printf("Total LHS Assembly Time: %f [s]\n",totLHSAssemblyTime);
}
