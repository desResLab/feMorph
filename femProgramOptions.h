#ifndef FEMPROGRAMOPTIONS_H
#define FEMPROGRAMOPTIONS_H

#include <string>

enum runModes{
  rmNORMAL,rmTRANSLATETOCVPRE,
  rmSIMPLEMAP,rmEXTRACTMESHQUALITY,
  rmMATCHFACELIST,rmMESHSKINTOCVPRE,
  rmCOMPUTEMODELEXPECTATIONS,
  rmCOMPUTEMODELWSS
};

class femProgramOptions{
  public:
    // CONSTRUCTOR
    femProgramOptions();
    // DATA MEMBERS
    runModes runMode;
    bool debugMode = false;
    bool reducedOutput = false;
    double tolerance;
    std::string inputFileName = "";
    std::string outputFileName = "";
    // MEMBER FUNCTIONS
    int getCommadLineOptions(int argc, char **argv);
};

#endif // FEMPROGRAMOPTIONS_H
