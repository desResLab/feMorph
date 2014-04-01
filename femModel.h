#ifndef FEMMODEL_H
#define FEMMODEL_H

#include <vector>
#include <string>

#include "femConstants.h"
#include "femNode.h"
#include "femElement.h"
#include "femFace.h"
#include "femEdge.h"
#include "femProperty.h"
#include "femModelSlice.h"
#include "femInputData.h"
#include "femResult.h"

class femModel
{
  public:
    std::vector<femNode*> nodeList;
    std::vector<femElement*> elementList;
    std::vector<femFace*> faceList;
    std::vector<femEdge*> edgeList;
    std::vector<femProperty*> propList;
    // Results
    std::vector<femResult*> resultList;
    // Enclosing Box and Model Centre
    double modelBox[6];
    double modelCentre[3];
    double weight = 0.0;
  public:
    // Constructor and Destructor
    femModel();
    ~femModel();
    // Other Member Functions
    void EvalModelBox();
    void EvalModelCentre();
    void FormElementFaceList();
    void FormElementEdgeList();
    void OrientateBoundaryFaceNodes();

    // Make Copies to other Model
    void CopyNodesTo(femModel* otherModel);
    void CopyElementsTo(femModel* otherModel);
    void CopyFacesTo(femModel* otherModel);
    void CopyPropertyTo(femModel* otherModel);

    // ====================
    // READ FUNCTIONALITIES
    // ====================
    // Read Whole Model From LS-DYNA File
    void ReadModelFromFile(std::string fileName);
    // Read Only Coords
    void ReadNodeCoordsFromFile(std::string fileName, bool skipFirstRow);
    // Read Only Element Connections
    void ReadElementConnectionsFromFile(std::string fileName, bool skipFirstRow);
    // Read Node Displacements From File
    void ReadNodeDisplacementsFromFile(std::string fileName, bool readRotations);
    // Read Node Coordinates From VTK Legacy
    void ReadModelNodesFromVTKFile(std::string fileName);
    // Read Element Connectivities From VTK Legacy
    void ReadModelElementsFromVTKFile(std::string fileName);
    // Read Model Results From VTK Legacy
    void ReadModelResultsFromVTKFile(std::string fileName);

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
    // Export to cvPRE
    void ExportToCvPre(double dispFactor, std::string pathName);
    void WriteCvPreFile(std::string fileName);
    void WriteAdjacenciesToFile(std::string fileName);
    void ExportBoundaryElementFiles(double dispFactor, int totalFaceGroups, std::string pathName);
    void ExportBoundaryNodeFiles(int totalFaceGroups, std::string pathName);
    void ExportElementFaceGroupToFile(std::string pathName, int groupID);
    void ExportNodeFaceGroupToFile(std::string pathName, int groupID);
    // Export Skin Faces of a given group to VTK   
    void ExportSkinFaceGroupToVTK(std::string fileName, double dispFactor, int groupNumber);
    // Write PolyFile for TETGEN
    void WriteSkinSMeshFile(std::string polyFileName);
    // Convert Node and Element File To CVPre
    int ConvertNodeAndElementsToCvPre(std::string nodeFileName, std::string elementFileName, bool skipFirstRow);

    // ==================================
    // CHECKS AND GEOMETRICAL EVALUATIONS
    // ==================================
    bool IsInsideLimits(double* nodeCoords);
    bool CheckNormalCompatibility(int firstElement, int secondElement);

    // =====================================
    // ASSIGMENT OF BOUNDARY FLOW CONDITIONS
    // =====================================
    void CreateBoundaryConditionFile(std::string inputFile);

    // ===============
    // MODEL ENQUIRIES
    // ===============
    double CheckMinimumElementVolume(double dispFactor);
    double CheckMinimumElementMixProduct(double dispFactor);
    void EvalModelQualityDistributions(std::string fileName, double* limitBox);

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
    // Group boundary faces based on normal
    void FormBoundaryFaceGroups(int &totalFaceGroups);
    // Form a model using the boundary faces only
    femModel* FormBoundaryFaceModel();
    // Assign 2D element property using normal
    void GroupFacesByNormal(int &currGroup);
    // Eval element Normal
    void evalElementNormal(int firstElement, double* normal);
    // Normalize Model Displacements
    void NormalizeDisplacements(double maxDisp);
    // Check If Two model have compatible Boxes
    bool isModelCompatible(femModel* other,double tolerance);

    // =======
    // MESHING
    // =======
    void MeshWithTetGen(std::string polyFileName);

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
    int FindEnclosingElementWithGrid(double dispFactor, double* nodeCoords, std::vector<int> &gridElementList);

    // ===================
    // RESULT MANIPULATION
    // ===================
    // Get Result Index From Label
    int getResultIndexFromLabel(std::string label);

    // ========================
    // Stenosis Parametrization
    // ========================
    // Seek the displacements factor giving a pre-set stonosis level
    double seekStenoticDisplacementFactor(femInputData* data, double targetStenosisLevel, bool debugMode);
    // Extract the stenosis for a given displacement factor
    double ExtractStenosisLevel(femInputData* data, double currDispFactor, std::vector<femModelSlice*> &slices, std::vector<double> &sliceAreas,
                                bool useDiameter, bool useOldDefinition, double* stenosisDef);
    // Slice Model Skin
    void SliceModelSkin(const int kStenosisSlices, double dispFactor, femInputData* data, std::vector<femModelSlice*> &slices);
    // GET AREA AT STENOSIS ORIGIN
    double GetReferenceArea(femInputData* data);
};

#endif // FEMMODEL_H
