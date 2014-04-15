#include "femProgramOptions.h"
#include "femUtils.h"

// CONSTRUCTOR
femProgramOptions::femProgramOptions(){
}

// ===============================================
// GET PROGRAM OPTIONS FROM COMMAND LINE ARGUMENTS
// ===============================================
int femProgramOptions::getCommadLineOptions(int argc, char **argv){

  // Declare
  int c;  
  FILE* echoFile;
  echoFile = fopen("options.echo","w");

  // Loop Through the Parameters
  while ((c = getopt (argc, argv, "f:o:ncmelsxbtdh")) != -1){
    switch (c){
      case 'f':
        inputFileName = std::string(optarg);
        fprintf(echoFile,"Input/Node File: %s\n",inputFileName.c_str());
        break;
      case 'o':
        outputFileName = std::string(optarg);
        fprintf(echoFile,"Output/Element File: %s\n",outputFileName.c_str());
        break;
      case 'n':
        runMode = rmNORMAL;
        fprintf(echoFile,"Run Mode: Normal\n");
        break;
      case 'c':
        runMode = rmTRANSLATETOCVPRE;
        fprintf(echoFile,"Run Mode: Translate to CVPre\n");
        break;
      case 'm':
        runMode = rmSIMPLEMAP;
        fprintf(echoFile,"Run Mode: Mapping\n");
        break;
      case 'e':
        runMode = rmEXTRACTMESHQUALITY;
        fprintf(echoFile,"Run Mode: Mesh Quality Extraction\n");
        break;
      case 'l':
        runMode = rmMATCHFACELIST;
        fprintf(echoFile,"Run Mode: Boundary face list matching\n");
        break;
      case 's':
        runMode = rmMESHSKINTOCVPRE;
        fprintf(echoFile,"Run Mode: Converting triangulation to CVPre Model\n");
        break;
      case 'x':
        runMode = rmCOMPUTEMODELEXPECTATIONS;
        fprintf(echoFile,"Run Mode: Computing Model Expectations\n");
        break;
      case 'b':
        runMode = rmCOMPUTEMODELWSS;
        fprintf(echoFile,"Run Mode: Computing Model WSS\n");
        break;
      case 't':
        tolerance = atof(optarg);
        fprintf(echoFile,"Tolerance: %s\n",optarg);
        break;
      case 'd':
        debugMode = true;
        fprintf(echoFile,"Running in Debug Mode\n");
        break;
      case 'h':
        femUtils::WriteAppHelp();
        break;
      case '?':
        if (optopt == 'f'){
          fprintf (echoFile, "Option -%c requires an input file name.\n", optopt);
        }else if (optopt == 'o'){
          fprintf (echoFile, "Option -%c requires an output file name.\n", optopt);
        }
      default:
        abort();
    }

    //for (int index = optind; index < argc; index++){
    //  fprintf(echoFile,"Non-option argument %s\n", argv[index]);
    //}
  }
  // Close Echo File
  fclose(echoFile);
  return 0;
}
