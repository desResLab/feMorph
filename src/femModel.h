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
#include "femProgramOptions.h"

using namespace std;

class femModel{
  public:
    vector<femNode*> nodeList;
    vector<femElement*> elementList;
    vector<femElement*> bcElementList;
    femDoubleVec bcElementValue;
    femDoubleMat bcElementNormal;
    vector<femFace*> faceList;
    vector<femEdge*> edgeList;
    vector<femProperty*> propList;
    // Parent element for each bc element
    femIntVec bcParentElement;
    // Element Velocity
    femDoubleMat elVelocity;
    // Source Array
    femIntVec sourceElement;
    femDoubleVec sourceValues;
    // Element Diffusivity
    femDoubleMat elDiffusivity;
    // Element Density
    femDoubleVec elDensity;
    // Element Viscosity
    femDoubleVec elViscosity;
    // Dirichelet BC Array
    femIntVec diricheletBCNode;
    femDoubleVec diricheletBCValues;
    // Neumann BC Array
    femIntVec neumannBCElement;
    femIntMat neumannBCFaceNodes;
    femDoubleVec neumannBCValues;
    // Maximum Number of Dofs per Node
    int maxNodeDofs;
    // Time integration parameters
    double timeStep;
    int totalSteps;
    int saveEvery;
    double alphaM;
    double alphaF;
    double gamma;
    // Solution Stages
    femIntVec solStages;
    // Use Prescribed Velocities
    bool usePrescribedVelocity;
    int prescribedVelType;
    // Storage for initial conditions
    femIntVec iniNodeNumbers;
    femIntVec iniDofNumber;
    femDoubleVec iniDofValue;

    // MPI Partitioning Information
    int totNodesInProc;
    int* localToGlobalNodes;

    // Results
    vector<femResult*> resultList;
    // Enclosing Box and Model Centre
    double modelBox[6];
    double modelCentre[3];
    double weight1;
    double weight2;

    // Constructor and Destructor
    femModel();
    virtual ~femModel();

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
    void ReadElementConnectionsFromFile(std::string fileName, bool skipFirstRow, bool numbersFromZero);
    // Read Node Displacements From File
    void ReadNodeDisplacementsFromFile(std::string fileName, bool readRotations);
    // Read File From VTK Legacy
    void ReadFromVTKLegacy(std::string fileName);
    // Read Node Coordinates From VTK Legacy
    void ReadModelNodesFromVTKFile(std::string fileName);
    // Read Element Connectivities From VTK Legacy
    void ReadModelElementsFromVTKFile(std::string fileName);
    // Read Model Results From VTK Legacy
    void ReadModelResultsFromVTKFile(std::string fileName);
    // Read Element Sources From File
    void ReadElementSourceFromFile(std::string fileName, bool skipFirstRow, bool numbersFromZero);
    // Read Dirichelet Boundary Conditions From File
    void ReadDirBCFromFile(std::string fileName, bool skipFirstRow, bool numbersFromZero);
    // Read Neumann Boundary Conditions From File
    void ReadNeumannBCFromFile(std::string fileName, bool skipFirstRow, bool numbersFromZero);
    // Read Element Diffusivity From File
    void ReadDiffusivityFromFile(std::string fileName, bool skipFirstRow, bool numbersFromZero);
    // Read Model From Text File
    void ReadFromFEMTextFile(std::string fileName);

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
    void ExportToCvPre(double dispFactor, std::string pathName, double angleLimit);
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
    int ConvertNodeAndElementsToCvPre(std::string nodeFileName, std::string elementFileName, bool vtkFile, bool skipFirstRow, double angleLimit);
    // Copy velocity Results to Vector
    void copyModelVelocitiesToVector(std::vector<std::vector<double> > &velocity);

    // ==================================
    // CHECKS AND GEOMETRICAL EVALUATIONS
    // ==================================
    bool IsInsideLimits(double* nodeCoords);
    bool CheckNormalCompatibility(int firstElement, int secondElement, double angleLimit);

    // =====================================
    // ASSIGMENT OF BOUNDARY FLOW CONDITIONS
    // =====================================
    void CreateBoundaryConditionFile(std::string inputFile);

    // ===============
    // MODEL ENQUIRIES
    // ===============
    double CheckMinimumElementVolume(double dispFactor);
    double CheckMinimumElementMixProduct(double dispFactor);
    void   EvalModelQualityDistributions(std::string fileName, double* limitBox);
    void   BuildParentElementList();

    // ==============
    // MODEL TOPOLOGY
    // ==============
    void getModelNodalTopology(femIntVec& diagPtr,femIntVec& rowPtr);

    // ==================
    // MODEL PARTITIONING
    // ==================
#ifdef USE_TRILINOS
    femModel* CreatePartition(int numPartitions);
#endif

    // ====================
    // MODEL MANIPULATUIONS
    // ====================
    // Fixed Element Connectivity
    void FixedElementConnectivities();
    // Perform Displacement Mapping
    void MapDisplacements(femProgramOptions* opts,femInputData* data,femModel* MappingModel,double dispScaleFactor);
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
    void FormBoundaryFaceGroups(int &totalFaceGroups, double angleLimit);
    // Form a model using the boundary faces only
    femModel* FormBoundaryFaceModel();
    // Assign 2D element property using normal
    void GroupFacesByNormal(int &currGroup, double angleLimit);
    // Eval element Normal: 2D Elements
    void eval2DElementNormal(int firstElement, double* normal);
    // Eval element Normal: 3D Elements
    void eval3DElementNormal(int elementID, int faceID, double* normal);
    // Normalize Model Displacements
    void NormalizeDisplacements(double maxDisp);
    // Check If Two model have compatible Boxes
    bool isModelCompatible(femModel* other,double tolerance);
    // Calculate Wall Shear Stresses
    void ComputeWSS();
    // Calculate Wall Shear Stresses Gradients
    void ComputeWSSGradients();
    // Apply Parametric Displacements
    void evalDisplacements(femInputData* data,double coordX,double coordY,double coordZ,double &dispX,double &dispY,double &dispZ);
    void ApplyParametricDisplacements(femInputData* data);

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
    // STENOSIS PARAMETRIZATION
    // ========================
    // Seek the displacements factor giving a pre-set stonosis level
    double seekStenoticDisplacementFactor(femInputData* data, double targetStenosisLevel, bool debugMode);
    // Extract the stenosis for a given displacement factor
    double ExtractStenosisLevel(femInputData* data, double currDispFactor, std::vector<femModelSlice*> &slices, std::vector<double> &sliceAreas,
                                bool useDiameter, bool useOldDefinition, double* stenosisDef);
    // Slice Model Skin
    void SliceModelSkin(const int kStenosisSlices, double dispFactor, femInputData* data, std::vector<femModelSlice*> &slices);
    // Get area at the origin of the stenosis
    double GetReferenceArea(femInputData* data);
};

#endif // FEMMODEL_H
