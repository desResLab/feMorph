#ifndef FEMMODELSEQUENCE_H
#define FEMMODELSEQUENCE_H

#include <string>
#include <vector>

#include "femModel.h"
#include "femTypes.h"

// SEQUENCE OF FEM MODELS
class femModelSequence{
  public:
    // CONSTRUCTOR
    femModelSequence();
    // DATA MEMBERS
    std::vector<femModel* > models;
    // MEMBER FUNCTIONS
    // IO
    void ReadFromWeightedListFile(std::string fileName);
    // COMPUTE WALL SHEAR STRESSES
    void ComputeWSS();
    // COMPUTE AV AND SD
    void ComputeResultStatistics(bool computeSD);
    void FixedElementConnectivities();
    void FormElementFaceList();
};

// Label Counter Struct
struct labelCounter{
  std::string label;
  int count;
  femResultType type;
};

#endif // FEMMODELSEQUENCE_H
