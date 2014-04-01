#include <string>

#include "femModelSequence.h"
#include "femUtils.h"

// CONSTRUCTOR
femModelSequence::femModelSequence(){

}

// READ MODEL SEQUENCE FROM WEIGHTED FILE LIST
void femModelSequence::ReadFromWeightedListFile(std::string fileName){
  // List with files and associated weights
  std::vector<femWeightedFileName*> fileList;
  std::string currFileName;
  double currWeight = 0.0;

  // Read files names and weights
  femUtils::ReadWeightedFileList(fileName,fileList);

  // Loop for all files in sequence
  for(unsigned int loopA=0;loopA<fileList.size();loopA++){
    currFileName = fileList[loopA]->fileName;
    currWeight = fileList[loopA]->weight;
    femModel* model = new femModel();
    // Read Nodes
    femUtils::WriteMessage("\n");
    femUtils::WriteMessage(std::string("Reading File ") + currFileName + "...\n");

    // Read nodes and elements only for the first file
    if(loopA == 0){
    femUtils::WriteMessage(std::string("Reading Nodes...\n"));
    model->ReadModelNodesFromVTKFile(currFileName);
    // Read Elements
    femUtils::WriteMessage(std::string("Reading Elements...\n"));
    model->ReadModelElementsFromVTKFile(currFileName);
    }

    // Read Results
    femUtils::WriteMessage(std::string("Reading Results...\n"));
    model->ReadModelResultsFromVTKFile(currFileName);
    // Assign Integration Weight
    model->weight = currWeight;
    // Add to models
    models.push_back(model);
  }
}

// FIND INDEX IN LABEL VECTOR IF LABEL AND TYPE ARE THE SAME
int isInLabelVector(std::string currLabel,femResultType currType,std::vector<labelCounter*> labelCount){
  bool found = false;
  int count = 0;
  while((!found)&&(count<labelCount.size())){
    found = ((labelCount[count]->label == currLabel)&&(labelCount[count]->type == currType));
    // Update Counter
    count++;
  }
  if(found){
    count--;
    return count;
  }else{
    return -1;
  }
}

// COMPUTE AV AND SD
void femModelSequence::ComputeResultStatistics(){

  // Count Result Labels in Models
  std::vector<labelCounter*> labelCount;
  std::string currLabel;
  femResultType currType;
  int currID = 0;
  int totEntities = 0;
  femModel* combModel = nullptr;
  int resIdx = 0;
  int currResult = 0;
  double currWeight = 0.0;

  for(size_t loopA=0;loopA<models.size();loopA++){
    for(size_t loopB=0;loopB<models[loopA]->resultList.size();loopB++){
      currLabel = models[loopA]->resultList[loopB]->label;
      currType = models[loopA]->resultList[loopB]->type;
      currID = isInLabelVector(currLabel,currType,labelCount);
      if(currID>-1){
        labelCount[currID]->count += 1;
      }else{
        // Add it to the labelCount
        labelCounter* lc = new labelCounter();
        lc->label = currLabel;
        lc->type = currType;
        lc->count = 1;
        labelCount.push_back(lc);
      }
    }
  }

  // If there is at least one label Create a new Model
  if(labelCount.size()>0){
    // Copy Nodes and Elements From First Model
    combModel = new femModel();
    models[0]->CopyNodesTo(combModel);
    models[0]->CopyElementsTo(combModel);
  }

  // Loop through all results that appear models.size() times
  for(int loopA=0;loopA<labelCount.size();loopA++){
    if(labelCount[loopA]->count == models.size()){

      // Combine all results with this Label
      // Average Values
      femResult* resAV = new femResult();
      resAV->label = labelCount[loopA]->label + "_AV";
      resAV->type = labelCount[loopA]->type;
      // Standard Deviations
      femResult* resSD = new femResult();
      resSD->label = labelCount[loopA]->label + "_SD";
      resSD->type = labelCount[loopA]->type;
      // Check Total Entities
      if(labelCount[loopA]->type == frNode){
        totEntities = combModel->nodeList.size();
      }else{
        totEntities = combModel->elementList.size();
      }
      // Loop through the entities
      for(int loopB=0;loopB<totEntities;loopB++){
        // Loop through all results to Compute the Average Value
        double avValue = 0.0;
        for(int loopC=0;loopC<models.size();loopC++){
          resIdx = models[loopC]->getResultIndexFromLabel(labelCount[loopA]->label);
          currResult = models[loopC]->resultList[resIdx]->values[loopB];
          currWeight = models[loopC]->weight;
          avValue += currResult*currWeight;
        }
        resAV->values.push_back(avValue);
        // Loop through all results to Compute the Standard Deviation
        double stdValue = 0.0;
        for(int loopC=0;loopC<models.size();loopC++){
          resIdx = models[loopC]->getResultIndexFromLabel(labelCount[loopA]->label);
          currResult = models[loopC]->resultList[resIdx]->values[loopB];
          currWeight = models[loopC]->weight;
          stdValue += (currResult-avValue)*(currResult-avValue)*currWeight;
        }
        resSD->values.push_back(sqrt(stdValue));
      }
      // Add both Results to The Model
      combModel->resultList.push_back(resAV);
      combModel->resultList.push_back(resSD);
    }
  }
  // Add Model To Sequence
  models.push_back(combModel);
}
