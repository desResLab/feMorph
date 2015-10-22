#ifndef FEMINTEGRATIONRULE_H
#define FEMINTEGRATIONRULE_H

#include "femTypes.h"

class femIntegrationRule{
  protected:
    intRuleType type;
  public:
    // DATA ACCESS
    int getTotGP(int totNodes,elDim dims);
    femDoubleMat getCoords(int totNodes,elDim dims);
    femDoubleVec getWeights(int totNodes,elDim dims);

    // CONSTRUCTOR
    femIntegrationRule(intRuleType type);
};

#endif // FEMINTEGRATIONRULE_H
