#ifndef FEMPROGRAMOPTIONS_H
#define FEMPROGRAMOPTIONS_H

#include <string>

using namespace std;

enum runModes{
  rmNORMAL,rmTRANSLATETOCVPRE,
  rmSIMPLEMAP,rmEXTRACTMESHQUALITY,
  rmMATCHFACELIST,rmMESHSKINTOCVPRE,
  rmCOMPUTEMODELEXPECTATIONS,
  rmCOMPUTEMODELWSS,rmSOLVEPOISSON,
  rmTESTELEMENTS,rmSOLVESTEADYSTATEADVECTIONDIFFUSION,
  rmSOLVEINCOMPRESSIBLENS,rmFFD
};

class femProgramOptions{
  public:
    // CONSTRUCTOR
    femProgramOptions();
    // DATA MEMBERS
    runModes runMode;
    bool debugMode;
    bool reducedOutput;
    // Use Legacy VTK as inputs
    bool useVTKFile;
    double tolerance;
    // Normal Angle for surface identification
    double angleLimit;
    string inputFileName;
    string outputFileName;

    // File Names
    string velocityFileName;
    string nodeFileName;
    string connectionFileName;
    string diffusivityFileName;
    string sourceFileName;
    string diricheletBCFileName;
    string neumannBCFileName;

    // Option for Advection Diffusion Solver
    int advDiffScheme;
    int advDiffVelType;
    int advDiffSourceType;

    // MEMBER FUNCTIONS
    int getCommadLineOptions(int argc, char **argv);
};

#endif // FEMPROGRAMOPTIONS_H
