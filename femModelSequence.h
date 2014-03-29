#ifndef FEMMODELSEQUENCE_H
#define FEMMODELSEQUENCE_H

class femModelSequence{
  public:
    // CONSTRUCTOR
    femModelSequence();
    // DATA MEMBERS
    std::vector<femModel*> models;
    // MEMBER FUNCTIONS
    // IO
    ReadFromWeightedListFile(std::string fileName);
    // COMPUTE AV AND SD
    ComputeStatistics();
};

#endif // FEMMODELSEQUENCE_H
