#include "femIntegrationRule.h"
#include "femException.h"

femIntegrationRule::femIntegrationRule(intRuleType type){
  std::vector<double> tempCoords;
  switch(type){
    case irFirstOrder:
      totalGP = 1;
      tempCoords.push_back(1.0/4.0);
      tempCoords.push_back(1.0/4.0);
      tempCoords.push_back(1.0/4.0);
      tempCoords.push_back(1.0/4.0);
      coords.push_back(tempCoords);
      weight.push_back(1.0/6.0);
      break;
    case irSecondOrder:
      throw femException("Not Implemented.\n");
      break;
  }
}
