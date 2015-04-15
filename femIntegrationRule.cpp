# include <math.h>

# include "femIntegrationRule.h"
# include "femException.h"
# include "femConstants.h"

femIntegrationRule::femIntegrationRule(intRuleType locType){
  type = locType;
}

// Get total number of Gauss Points
int femIntegrationRule::getTotGP(int totNodes){
  int result;
  switch(type){
    case irFirstOrder:
      result = 1;
      break;
    case irSecondOrder:
      if(totNodes == kRodNodes){
        result = 2;
      }else if(totNodes == kTri3Nodes){
        result = 3;
      }else if(totNodes == kTetra4Nodes){
        result = 4;
      }else if(totNodes == kTetra10Nodes){
        result = 4;
      }else if(totNodes == kHexa8Nodes){
        result = 8;
      }else{
        throw femException("ERROR: Integration Rule not supported!");
      }
      break;
  }
  // Return Res
  return result;
}

// Get Intergration Rule Coordinates
femDoubleMat femIntegrationRule::getCoords(int totNodes){
  femDoubleVec tempCoords;
  femDoubleMat result;
  switch(type){
    // First Order Integration Rules
    case irFirstOrder:
      if(totNodes == kTri3Nodes){
        // First Order Integration for Triangles
        tempCoords.push_back(1.0/3.0);
        tempCoords.push_back(1.0/3.0);
        tempCoords.push_back(1.0/3.0);
        result.push_back(tempCoords);
      }else if(totNodes == kTetra4Nodes){
        // First Order Integration for Tetrahedrals
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        result.push_back(tempCoords);
      }else if(totNodes == kTetra10Nodes){
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        tempCoords.push_back(1.0/4.0);
        result.push_back(tempCoords);
      }else if(totNodes == kHexa8Nodes){
        tempCoords.push_back(0.0);
        tempCoords.push_back(0.0);
        tempCoords.push_back(0.0);
        result.push_back(tempCoords);
      }else{
        throw femException("ERROR: Integration Rule not supported!");
      }
      break;
    // Second Order Integration Rules
    case irSecondOrder:
      if(totNodes == kTri3Nodes){
        // Point 1
        tempCoords.clear();
        tempCoords.push_back(0.5);
        tempCoords.push_back(0.5);
        tempCoords.push_back(0.0);
        result.push_back(tempCoords);
        // Point 2
        tempCoords.clear();
        tempCoords.push_back(0.0);
        tempCoords.push_back(0.5);
        tempCoords.push_back(0.5);
        result.push_back(tempCoords);
        // Point 3
        tempCoords.clear();
        tempCoords.push_back(0.5);
        tempCoords.push_back(0.0);
        tempCoords.push_back(0.5);
        result.push_back(tempCoords);
      }else if(totNodes == kTetra4Nodes){
        // Point 1
        tempCoords.clear();
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 2
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 3
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 4
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        result.push_back(tempCoords);
      }else if(totNodes == kTetra10Nodes){
        // Point 1
        tempCoords.clear();
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 2
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 3
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        tempCoords.push_back(0.13819660);
        result.push_back(tempCoords);
        // Point 4
        tempCoords.clear();
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.13819660);
        tempCoords.push_back(0.58541020);
        result.push_back(tempCoords);
      }else if(totNodes == kHexa8Nodes){
        // Point 1
        tempCoords.clear();
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 2
        tempCoords.clear();
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 3
        tempCoords.clear();
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 4
        tempCoords.clear();
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 5
        tempCoords.clear();
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 6
        tempCoords.clear();
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 7
        tempCoords.clear();
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        result.push_back(tempCoords);
        // Point 8
        tempCoords.clear();
        tempCoords.push_back(-1.0/sqrt(3.0));
        tempCoords.push_back(1.0/sqrt(3.0));
        tempCoords.push_back(-1.0/sqrt(3.0));
        result.push_back(tempCoords);
      }else{
        throw femException("ERROR: Integration Rule not supported!");
      }
    }
  // Return Res
  return result;
}

// Get Integration Rule Weights
femDoubleVec femIntegrationRule::getWeights(int totNodes){
  femDoubleVec result;
  double origVolume = 0.0;
  switch(type){
    case irFirstOrder:
      if(totNodes == kTri3Nodes){
        // Second Order Integration Rule for Triangles
        origVolume = 0.5;
        result.push_back(origVolume*1.0);
      }else if(totNodes == kTetra4Nodes){
        origVolume = 1.0/6.0;
        result.push_back(origVolume*1.0);
      }else if(totNodes == kTetra10Nodes){
        origVolume = 1.0/6.0;
        result.push_back(origVolume*1.0);
      }else if(totNodes == kHexa8Nodes){
        origVolume = 8.0;
        result.push_back(origVolume*1.0);
      }else{
        throw femException("ERROR: Integration Rule not supported!");
      }
      break;
    case irSecondOrder:
      if(totNodes == kTri3Nodes){
        // Second Order Integration Rule for Triangles
        origVolume = 0.5;
        result.push_back(origVolume*(1.0/3.0));
        result.push_back(origVolume*(1.0/3.0));
        result.push_back(origVolume*(1.0/3.0));
      }else if(totNodes == kTetra4Nodes){
        origVolume = 1.0/6.0;
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
      }else if(totNodes == kTetra10Nodes){
        origVolume = 1.0/6.0;
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
        result.push_back(origVolume*(1.0/4.0));
      }else if(totNodes == kHexa8Nodes){
        origVolume = 8.0;
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
        result.push_back(origVolume*(1.0/8.0));
      }else{
        throw femException("ERROR: Integration Rule not supported!");
      }
      break;
    }
    // Return Res
    return result;
}


