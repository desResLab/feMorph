#ifndef FEMMODEL_H
#define FEMMODEL_H

#include <vector>
#include <string>

#include "femNode.h"
#include "femElement.h"
#include "femFace.h"
#include "femProperty.h"
#include "femInputData.h"

class femModel
{
  public:
    std::vector<femNode*> nodeList;
    std::vector<femElement*> elementList;
    std::vector<femFace*> faceList;
    std::vector<femProperty*> propList;
    double* modelBox;
  public:
    // Constructor and Destructor
    femModel();
    ~femModel();
    // Other Member Functions
    void EvalModelBox();
    void FormElementFaceList();

    // Make Copies to other Model
    void CopyElementsTo(femModel* otherModel);
    void CopyFacesTo(femModel* otherModel);
    void CopyPropertyTo(femModel* otherModel);

    // ==========
    // NODE TOOLS
    // ==========
    bool IsOutsideLimits(double* nodeCoords);

    // =============
    // ELEMENT TOOLS
    // =============
    void evalElementCentroid(int elementID, double* centroid);
    double evalPointToElementDistance(int elementID, double* pointCoords);

    // ====================
    // READ FUNCTIONALITIES
    // ====================
    // Read Whole Model From File
    void ReadModelFromFile(std::string fileName);
    // Read Only Coords
    void ReadNodeCoordsFromFile(std::string fileName);
    // Read Only Element Connections
    void ReadElementConnectionsFromFile(std::string fileName);
    // Read Node Displacements From File
    void ReadNodeDisplacementsFromFile(std::string fileName, bool readRotations);

    // =====================
    // WRITE FUNCTIONALITIES
    // =====================
    // Write Coords To File
    void WriteNodeCoordsToFile(std::string fileName);
    // Write Element Connections To File
    void WriteElementConnectionsToFile(std::string fileName);
    // Node List Manipulation
    int GetNodeIDFromNumber(int number);

    // ====================
    // MODEL MANIPULATUIONS
    // ====================
    // Find Element Containing a Given Node
    int FindEnclosingElement(double* nodeCoords);
    // Interpolate Element Displacement Field at a given Node
    void InterpolateElementDisplacements(double* nodeCoords, int elementID, double* nodeDisps);
    // Perform Displacement Mapping
    void MapDisplacements(femModel* MappingModel, femInputData* data, double dispScaleFactor);
    // Get next element when finding the enclosing one
    void getNextElement(int currElement, double* nodeCoords, int &nextElement, double &nextDistance);
    // Transform Model by translating/scaling/rotations
    femModel* TransformModel(femInputData* data, double* stenosisBox);
    // Get Stenosis Box
    void GetStenosisBox(femInputData* data, double* limRect);
    // Rotate Model
    void RotateModel(double angle1, double* axis1);
};

#endif // FEMMODEL_H
