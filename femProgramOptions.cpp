#include "femProgramOptions.h"
#include "femUtils.h"

femProgramOptions::femProgramOptions(){
}

// ===============================================
// GET PROGRAM OPTIONS FROM COMMAND LINE ARGUMENTS
// ===============================================
int femProgramOptions::getCommadLineOptions(int argc, char **argv){
  int c;
  while ((c = getopt (argc, argv, "f:o:ndcme")) != -1){
    switch (c){
      case 'f':
        inputFileName = optarg;
        break;
      case 'o':
        outputFileName = optarg;
        break;
      case 'n':
        runMode = rmNORMAL;
        break;
      case 't':
        runMode = rmTRANSLATETOCVPRE;
        break;
      case 'm':
        runMode = rmSIMPLEMAP;
        break;
      case 'e':
        runMode = rmEXTRACTMESHQUALITY;
        break;
      case 'l':
        runMode = rmMATCHFACELIST;
        break;
      case 's':
        runMode = rmMESHSKINTOCVPRE;
        break;
      case 'x':
        runMode = rmCOMPUTEMODELEXPECTATIONS;
        break;
      case 'd':
        debugMode = true;
        break;
      case 'h':
        femUtils::WriteAppHelp();
        break;
      case '?':
        if (optopt == 'f'){
          fprintf (stderr, "Option -%c requires an input file name.\n", optopt);
        }else if (optopt == 'o'){
          fprintf (stderr, "Option -%c requires an output file name.\n", optopt);
        }
      default:
        abort ();
    }

    for (int index = optind; index < argc; index++){
      printf ("Non-option argument %s\n", argv[index]);
    }
  }
  return 0;
}
