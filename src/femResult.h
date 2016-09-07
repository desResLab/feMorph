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
    int numComponents;
    femDoubleMat values;
};

#endif // FEMRESULT_H
