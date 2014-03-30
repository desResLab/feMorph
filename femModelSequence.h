#ifndef FEMMODELSEQUENCE_H
#define FEMMODELSEQUENCE_H

#include <string>
#include <vector>

#include "femModel.h"
#include "femTypes.h"


// Label Counter Struct
struct labelCounter{
  std::string label = "";
  int count = 0;
  femResultType type;
};

// SEQUENCE OF FEM MODELS
class femModelSequence{
  public:
    // CONSTRUCTOR
    femModelSequence();
    // DATA MEMBERS
    std::vector<femModel*> models;
    // MEMBER FUNCTIONS
    // IO
    void ReadFromWeightedListFile(std::string fileName);
    // COMPUTE AV AND SD
    void ComputeResultStatistics();
};

#endif // FEMMODELSEQUENCE_H
