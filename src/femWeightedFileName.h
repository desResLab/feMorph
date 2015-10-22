#ifndef FEMWEIGHTEDFILENAME_H
#define FEMWEIGHTEDFILENAME_H

#include <string>

class femWeightedFileName{
  public:
    // CONSTRUCTOR
    femWeightedFileName(std::string currfileName, double currWeight1, double currWeight2);
    // DATA MEMBERS
    std::string fileName;
    double weight1;
    double weight2;
};

#endif // FEMWEIGHTEDFILENAME_H
