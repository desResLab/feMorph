#ifndef FEMMODEL_H
#define FEMMODEL_H

class femModel
{
  protected:
    vector<femNode> nodeList;
    vector<femElement> elementList;
    vector<femProperty> propList;
  
  public:
    // Constructor and Destructor
    femModel();
    ~femModel();
    // Other Member Functions
    void ReadModelFromLSDYNAFile();
};

#endif // FEMMODEL_H
