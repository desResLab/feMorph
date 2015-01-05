#ifndef FEMPROGRAMOPTIONS_H
#define FEMPROGRAMOPTIONS_H

#include <string>

using namespace std;

enum runModes{
  rmNORMAL,rmTRANSLATETOCVPRE,
  rmSIMPLEMAP,rmEXTRACTMESHQUALITY,
  rmMATCHFACELIST,rmMESHSKINTOCVPRE,
  rmCOMPUTEMODELEXPECTATIONS,
  rmCOMPUTEMODELWSS,rmSOLVEPOISSON
};

class femProgramOptions{
  public:
    // CONSTRUCTOR
    femProgramOptions();
    // DATA MEMBERS
    runModes runMode;
    bool debugMode = false;
    bool reducedOutput = false;
    // Use Legacy VTK as inputs
    bool useVTKFile = false;
    double tolerance;
    // Normal Angle for surface identification
    double angleLimit;
    string inputFileName = "";
    string outputFileName = "";

    // File Names
    string nodeFileName = "";
    string connectionFileName = "";
    string sourceFileName = "";
    string diricheletBCFileName = "";
    string neumannBCFileName = "";

    // MEMBER FUNCTIONS
    int getCommadLineOptions(int argc, char **argv);
};

#endif // FEMPROGRAMOPTIONS_H
