#ifndef FEMINTEGRATIONRULE_H
#define FEMINTEGRATIONRULE_H

#include "femTypes.h"

class femIntegrationRule{
  public:
    int totalGP;
    femDoubleMat coords;
    femDoubleVec weight;
    // CONSTRUCTOR
    femIntegrationRule(intRuleType type);
};

#endif // FEMINTEGRATIONRULE_H
