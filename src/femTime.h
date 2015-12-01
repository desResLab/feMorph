#ifndef FEMTIME_H
#define FEMTIME_H

# include <stdio.h>

class femTime{
  public:
    // Data Members
    static double LocalSFTime;
    static double JacobianMatTime;
    static double JacInversionTime;
    static double globShDerivAllocTime;
    static double totLHSAssemblyTime;

    // Member Functions
    static void printToScreen();
};

#endif // FEMTIME_H
