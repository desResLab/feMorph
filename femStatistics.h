#ifndef FEMSTATISTICS_H
#define FEMSTATISTICS_H

#include "femModel.h"

namespace femStatistics
{

// Normalize Bin Array
inline void NormalizeBinArray(int size, double* binArray,double currInterval){
  double sum = 0.0;
  // Compute Summation
  for(int loopA=0;loopA<size;loopA++){
    sum += binArray[loopA];
  }
  // Exit if zero sum
  if (fabs(sum)<kMathZero){
    return;
  }
  // Normalize
  for(int loopA=0;loopA<size;loopA++){
    binArray[loopA] /= (sum*currInterval);
  }
}

// Assign to BIN
inline void AssignToBin(double currValue, int numberOfBins, double* binMin, double* binMax, double* binArray){
  bool found = false;
  int count = 0;
  bool isMoreThanMin = false;
  bool isLessThanMax = false;
  while ((!found)&&(count<numberOfBins)){
    if (fabs(currValue-binMin[0])<kMathZero){
      isMoreThanMin = (currValue >= binMin[count] - kMathZero);
    }else{
      isMoreThanMin = (currValue > binMin[count]);
    }
    if (fabs(currValue-binMin[numberOfBins-1])<kMathZero){
      isLessThanMax = (currValue <= binMax[count]);
    }else{
      isLessThanMax = (currValue <= binMax[count] + kMathZero);
    }
    found = (isMoreThanMin)&&(isLessThanMax);
    // Update
    if (!found){
      count++;
    }
  }
  if (found){
    // Increase Bin Count
    binArray[count] = binArray[count] + 1.0;
  }else{
    throw femException("Error: Value Cannot fit in Bin.\n");
  }
}

// FORM BIN LIMITS FOR SINGLE SCAN
inline void FormBinLimits(femModel* model, int modelQuality, double &currInterval, int numberOfBins, double* binMin, double* binMax, double* binCenter, double* limitBox){
  // Initialize Limits
  double  minRange = std::numeric_limits<double>::max();
  double  maxRange = -std::numeric_limits<double>::max();
  double  currValue = 0.0;
  double centroid[3] = {0.0};
  for(unsigned int loopA=0;loopA<model->elementList.size();loopA++){
    // Get Element Centroid
    model->elementList[loopA]->evalElementCentroid(model->nodeList,centroid);
    // Check if inside Limits
    if(femUtils::isInsideLimits(centroid,limitBox)){
      // Get quantity
      if(modelQuality == kVolume){
        currValue = model->elementList[loopA]->EvalVolume(0.0,model->nodeList);
      }else if(modelQuality == kMixedProduct){
        currValue = model->elementList[loopA]->EvalMixProduct(model->nodeList);
      }
      // Assign Values
      if(currValue>maxRange){
        maxRange = currValue;
      }
      if(currValue<minRange){
        minRange = currValue;
      }
    }
  }
  // If minRange and maxRange are the same than add something
  if (fabs(maxRange - minRange)<kMathZero){
    minRange = minRange - 1.0;
    maxRange = maxRange + 1.0;
  }
  // Fill the bin arrays
  double currPtr = minRange;
  currInterval = ((maxRange - minRange)/(double)numberOfBins);
  for(int loopA=0;loopA<numberOfBins;loopA++){
    binMin[loopA] = currPtr;
    binMax[loopA] = currPtr + currInterval;
    binCenter[loopA] = 0.5*(binMin[loopA] + binMax[loopA]);
    // Update
    currPtr += currInterval;
  }
}

} // femStatistics


#endif // FEMSTATISTICS_H
