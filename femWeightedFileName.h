#ifndef FEMWEIGHTEDFILENAME_H
#define FEMWEIGHTEDFILENAME_H

#include <string>

class femWeightedFileName
{
public:
    // CONSTRUCTOR
    femWeightedFileName(std::string currfileName, double currWeight);
    // DATA MEMBERS
    std::string fileName;
    double weight;
};

#endif // FEMWEIGHTEDFILENAME_H
