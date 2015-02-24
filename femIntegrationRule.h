#ifndef FEMINTEGRATIONRULE_H
#define FEMINTEGRATIONRULE_H

#include "femTypes.h"

class femIntegrationRule{
  protected:
    intRuleType type;
  public:
    // DATA ACCESS
    int getTotGP(int totNodes);
    femDoubleMat getCoords(int totNodes);
    femDoubleVec getWeights(int totNodes);

    // CONSTRUCTOR
    femIntegrationRule(intRuleType type);
};

#endif // FEMINTEGRATIONRULE_H
