#ifndef FEMRESULT_H
#define FEMRESULT_H

#include <string>
#include <vector>

#include "femTypes.h"

// FEM RESULT - GENERIC
class femResult{
  public:
    femResult();
    // DATA MEMBER
    std::string label;
    femResultType type;
    std::vector<double> values;
};

#endif // FEMRESULT_H
