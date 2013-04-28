#ifndef FEMMODEL_H
#define FEMMODEL_H

#include <vector>
#include <string>

#include "femConstants.h"
#include "femNode.h"
#include "femElement.h"
#include "femFace.h"
#include "femProperty.h"
#include "femModelSlice.h"
#include "femInputData.h"

class femModel
{
  public:
    std::vector<femNode*> nodeList;
    std::vector<femElement*> elementList;
    std::vector<femFace*> faceList;
    std::vector<femFace*> boundaryFaceList;
    std::vector<femProperty*> propList;
    double modelBox[6];
    double modelCentre[3];
  public:
    // Constructor and Destructor
    femModel();
    ~femModel();
    // Other Member Functions
    void EvalModelBox();
    void EvalModelCentre();
    void FormElementFaceList();
    void FormBoundaryFaceList();

    // Make Copies to other Model
    void CopyElementsTo(femModel* otherModel);
    void CopyFacesTo(femModel* otherModel);
    void CopyPropertyTo(femModel* otherModel);

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
    void WriteNodeCoordsToFile(double dispFactor, std::string fileName);
    // Write Element Connections To File
    void WriteElementConnectionsToFile(std::string fileName);
    // Node List Manipulation
    int GetNodeIDFromNumber(int number);
    // Export Model to VTK Legacy
    void ExportToVTKLegacy(std::string fileName);

    // ======
    // CHECKS
    // ======
    bool IsInsideLimits(double* nodeCoords);

    // ====================
    // MODEL MANIPULATUIONS
    // ====================
    // Perform Displacement Mapping
    void MapDisplacements(femModel* MappingModel, femInputData* data, double dispScaleFactor);
    // Get next element when finding the enclosing one
    void getNextElement(int currElement, double* nodeCoords, int &nextElement, double &nextDistance);
    // Get Stenosis Box
    void GetStenosisBox(femInputData* data, double* limRect);
    // Transform Model by translating/scaling/rotations
    femModel* TransformModel(femInputData* data, double* stenosisBoxCenter, double* stenosisBox);
    // Rotate Model
    void RotateModel(double angle, double* axis);
    // Write Rotate Node List from Box
    void createStenosisNodeList(double* stenosisBox, femInputData* data, std::vector<femNode*> &steNodeList);

    // ========================
    // ENCLOSING ELEMENT SEARCH
    // ========================
    // Find Element Containing a Given Node: Presearch
    int FindEnclosingElementPre(double* nodeCoords);
    // Find Element Containing a Given Node: Walking Algortihm
    int FindEnclosingElementPost(double* nodeCoords, int startingElement);
    // Find Element Containing a Given Node: Use Adjacencies
    int FindEnclosingElementWithAdj(double* nodeCoords);
    // Find Element Containing a Given Node: Use Grid
    int FindEnclosingElementWithGrid(double* nodeCoords, std::vector<int> &gridElementList);

    // ========================
    // Stenosis Parametrization
    // ========================
    // Seek the displacements factor giving a pre-set stonosis level
    double seekStenoticDisplacementFactor(femInputData* data, double targetStenosisLevel, bool debugMode);
    // Extract the stenosis for a given displacement factor
    double ExtractStenosisLevel(femInputData* data, double currDispFactor, std::vector<femModelSlice*> &slices, std::vector<double> &sliceAreas);
    // Slice Model Skin
    void SliceModelSkin(const int kStenosisSlices, double dispFactor, femInputData* data, std::vector<femModelSlice*> &slices);
};

#endif // FEMMODEL_H
