#ifndef FEMFFDOPTIONS_H
#define FEMFFDOPTIONS_H

#include <iostream>
#include <fstream>
#include <string>

# include "femTypes.h"
# include <boost/algorithm/string.hpp>

struct ffdGrid{
  double minPoint[3];
  int totPoints[3];
  double gridAxis_S[3];
  double gridAxis_T[3];
  double gridAxis_U[3];
  femIntVec dispNodes;
  femDoubleMat dispVals;
};

class femFFDOptions{
public:
  // Data Members
  string inputFileName;
  std::vector<ffdGrid> ffdData;
  // Member Functions
  femFFDOptions();
  ~femFFDOptions();
  // Read the FFD input file
  void readFromFile(string fileName);
};

#endif // FEMFFDOPTIONS_H
