#include "femElement.h"
#include "femNode.h"

#include "femException.h"
#include "femConstants.h"
#include "femUtils.h"

// Constructor
femElement::femElement(int number, int prop, int totalNodes, int* connections)
{
  // Assign Initial Values
  elementNumber = number;
  propertyNumber = prop;
  for(int loopA=0;loopA<totalNodes;loopA++){
    elementConnections.push_back(connections[loopA]);
  }
}

// Copy Constructor
femElement::femElement(femElement* other){
    // Numbers
    elementNumber = other->elementNumber;
    propertyNumber = other->propertyNumber;
    // Element Connectioes
    for(unsigned int loopA=0;loopA<other->elementConnections.size();loopA++){
      elementConnections.push_back(other->elementConnections[loopA]);
    }
    // Element Faces
    for(unsigned int loopA=0;loopA<other->elementFaces.size();loopA++){
      elementFaces.push_back(other->elementFaces[loopA]);
    }
}

femElement::~femElement()
{
}

// ===============================
// Eval TETRA10 Volume Coordinates
// ===============================
void femTetra10::EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){
  int* connections= new int[4];

  // Copy the first four connctions !!! Complete...
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    connections[loopA] = elementConnections[loopA];
  }

  // Create a temporay Tet4 Element with the First 4 nodes
  femTetra4* tet4 = new femTetra4(1,1,kTetra4Nodes,connections);
  // Compute Final Volume Coordinates
  double tet4VolCoords[kTetra4Nodes] = {0.0};
  tet4->EvalVolumeCoordinates(pointCoords,nodeList,tet4VolCoords);

  // Compute The Area Coordinates for the Full Quadratic Tetrahedron
  volCoords[0] = tet4VolCoords[0]*(2.0*tet4VolCoords[0]-1.0);
  volCoords[1] = tet4VolCoords[1]*(2.0*tet4VolCoords[1]-1.0);
  volCoords[2] = tet4VolCoords[2]*(2.0*tet4VolCoords[2]-1.0);
  volCoords[3] = tet4VolCoords[3]*(2.0*tet4VolCoords[3]-1.0);
  volCoords[4] = 4.0*tet4VolCoords[0]*tet4VolCoords[1];
  volCoords[5] = 4.0*tet4VolCoords[1]*tet4VolCoords[2];
  volCoords[6] = 4.0*tet4VolCoords[2]*tet4VolCoords[0];
  volCoords[7] = 4.0*tet4VolCoords[0]*tet4VolCoords[3];
  volCoords[8] = 4.0*tet4VolCoords[1]*tet4VolCoords[3];
  volCoords[9] = 4.0*tet4VolCoords[2]*tet4VolCoords[3];

  // Delete Temporary Tet4 element
  delete tet4;
}

// =======================================
// Assemble Tetra4 Coordinates in a Matrix
// =======================================
void femTetra4::AssembleTetCoordsMat(std::vector<femNode*> &nodeList, double** coordMat){
  // Fill Matrix
  int currNode = 0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    // Get Current Node
    currNode = elementConnections[loopA];
    // Store In Matrix
    for(int loopB=0;loopB<3;loopB++){
      coordMat[loopB][loopA] = nodeList[currNode]->coords[loopB];
    }
  }
}

// ====================================
// Eval Tetra4 Volume from Coord Matrix
// ====================================
double EvalTetVolumeFromCoordMat(double** coordMat){
  // Allocate
  double firstVec[3] = {0.0};
  double secondVec[3] = {0.0};
  double thirdVec[3] = {0.0};

  // Get The Vectors
  // First
  firstVec[0] = coordMat[0][1] - coordMat[0][0];
  firstVec[1] = coordMat[1][1] - coordMat[1][0];
  firstVec[2] = coordMat[2][1] - coordMat[2][0];
  // Second
  secondVec[0] = coordMat[0][2] - coordMat[0][0];
  secondVec[1] = coordMat[1][2] - coordMat[1][0];
  secondVec[2] = coordMat[2][2] - coordMat[2][0];
  // Third
  thirdVec[0] = coordMat[0][3] - coordMat[0][0];
  thirdVec[1] = coordMat[1][3] - coordMat[1][0];
  thirdVec[2] = coordMat[2][3] - coordMat[2][0];

  // Eval External Product
  double auxVec[3] = {0.0};
  femUtils::Do3DExternalProduct(secondVec,thirdVec,auxVec);

  // Eval Internal Product
  return (femUtils::Do3DInternalProduct(firstVec,auxVec)/6.0);
}

// ========================
// Get External Tet Volumes
// ========================
void EvalExternalTetVolumes(double* pointCoords, double** coordMat, double* extTetVol){
  // Allocate current coordinate matrix
  double** currCoordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    currCoordMat[loopA] = new double[kTetra4Nodes];
  }
  // Eval Volume Coords
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    for(int loopB=0;loopB<kTetra4Nodes;loopA++){
      if (loopB != loopA){
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = coordMat[loopC][loopB];
        }
      }else{
        for(int loopC=0;loopC<3;loopC++){
          // Assemble Modified Coord Matrix With adHoc Coords
          currCoordMat[loopC][loopB] = pointCoords[loopC];
        }
      }
    }
    // Eval Volume Using CurrCoordMat
    extTetVol[loopA] = EvalTetVolumeFromCoordMat(currCoordMat);
  }
  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] currCoordMat[loopA];
  }
  delete [] currCoordMat;
}



// ========================================
// Eval Volume Coordinates of Point in Tet4
// ========================================
void femTetra4::EvalVolumeCoordinates(double* pointCoords, std::vector<femNode*> &nodeList, double* volCoords){

  // Put The Coordinates in a Matrix
    // Allocate current coordinate matrix
  double** coordMat = new double*[3];
  for(int loopA=0;loopA<3;loopA++){
    coordMat[loopA] = new double[kTetra4Nodes];
  }
  AssembleTetCoordsMat(nodeList,coordMat);

  // Eval Tethrahedral Volume
  double TetVolume = EvalTetVolumeFromCoordMat(coordMat);

  // Eval Other Volumes
  double currVolumes[kTetra4Nodes] = {0.0};
  EvalExternalTetVolumes(pointCoords,coordMat,currVolumes);

  // Eval Shape Function Summation
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++)
  {
    sum += currVolumes[loopA];
  }

  // Compute Final Volume Coordinates
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    volCoords[loopA] = (currVolumes[loopA]/TetVolume);
  }

  // Deallocate
  for(int loopA=0;loopA<3;loopA++){
    delete [] coordMat[loopA];
  }
  delete [] coordMat;
}


// Check if Node is Inside
bool femTetra4::isNodeInsideElement(double* pointCoords,std::vector<femNode*> &nodeList){
  // Init Result
  bool isInside = false;

  // Eval Volume Coordinates
  double volCoords[kTetra4Nodes] = {0.0};
  EvalVolumeCoordinates(pointCoords,nodeList,volCoords);

  // Compute Result
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    isInside = (isInside)&&(volCoords[loopA] >= 0.0);
  }

  // Check if the volume coordinates sum up to one
  double sum = 0.0;
  for(int loopA=0;loopA<kTetra4Nodes;loopA++){
    sum += volCoords[loopA];
  }
  if (fabs(sum-1.0)>kMathZero){
    throw new femException("Error: Internal. Tet4 Shape Function don't sum up to one");
  }
  // Return
  return isInside;
}

// =================================
// Interpolate Element Displacements
// =================================
void femElement::InterpolateElementDisplacements(double* nodeCoords, std::vector<femNode*> &nodeList, double* intDisps){

  double* volCoords = new double[int(elementConnections.size())];
  // Eval Shape Function Values for Current Node
  EvalVolumeCoordinates(nodeCoords,nodeList,volCoords);

  // Loop through the nodes
  for(int loopA=0;loopA<kNodeDofs;loopA++){
    intDisps[loopA] = 0.0;
  }

  // Loop through the nodes
  int currNode = 0;
  double* currNodeDisps = nullptr;
  for(unsigned int loopA=0;loopA<elementConnections.size();loopA++){
    // Get Node Number
    currNode = elementConnections[loopA];
    // Get Node Displacements
    currNodeDisps = nodeList[currNode]->getNodeDisplacements();
    // Loop through the node dofs
    for(int loopB=0;loopB<kNodeDofs;loopB++){
      intDisps[loopB] += currNodeDisps[loopB]*volCoords[loopA];
    }
  }
}


/*// Transform Model
Function TransformNode(MappingModelID: Integer4;
                       XYZ: Array3Doubles;
                       ParamValue: Double;
                       Var NewXYZ: Array3Doubles): Integer4;
Const
  CurrentResultCase = 1;
Var
  LoopA,LoopB: Integer4;
  Error: Integer4;
  Count: Integer4;
  Found: Boolean;
  CoordMat: Double2DArray;
  FoundElement: Integer4;
  Connection: ConnectionArray;
  MorphedDisp: Array3Doubles;
  NodeRes: NodeResultArray;
  TotalMappingElements: Integer4;
  VolCoords: Array10Doubles;
  Sum: Double;
  VCoords: Array4Doubles;
Begin
  // Init Result
  Result:=0;

  // Get The total Number Of Mapping elements
  Error:=St7GetTotal(MappingModelID,tyBRICK,TotalMappingElements);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Find the Corresponding Element: Full Search: TEMPORARY
  Found:=FALSE;
  Count:=0;
  While (Not(Found))And(Count<TotalMappingElements) Do
  Begin
    // Update Counter
    Inc(Count);
    // Check If foundEvalExternalTetVolumes
    Found:=IsPointInsideLinearTet(MappingModelID,XYZ,Count);
    If Found Then
  End;

  // Check If found
  If (Found) Then
  Begin
    // Store the Element where the node belongs to
    FoundElement:=Count;

    // Eval Volume Coordinates for Quandratic Tetra
    VolCoords:=EvalTet10VolCoordinates(MappingModelID,XYZ,FoundElement);

    // Get The element Connections for the Node
    Error:=St7GetElementConnection(MappingModelID,tyBRICK,FoundElement,Connection);
    If (Error<>ERR7_NoError) Then FatalError(Error);

    // Initialize Results
    For LoopA:=0 To 2 Do MorphedDisp[LoopA]:=0.0;

    // Interpolate the displacement Field
    For LoopA:=0 To (Connection[0]-1) Do
    Begin
      // Get Node Displacement For Transformed Element
      Error:=St7GetNodeResult(MappingModelID,kNodeDisp,Connection[LoopA+1],CurrentResultCase,NodeRes);
      If (Error<>ERR7_NoError) Then FatalError(Error);

      // Interpolate Morphed Displacements
      For LoopB:=0 To 2 Do MorphedDisp[LoopB]:=MorphedDisp[LoopB]+VolCoords[LoopA]*NodeRes[LoopB]*ParamValue;
    End;

    // Store the New Coordinates
    For LoopA:=0 To 2 Do NewXYZ[LoopA]:=XYZ[LoopA]+MorphedDisp[LoopA];
  End Else Begin

    // If Cannot Find the Element than Use the Same Coordinates and Exit
    For LoopA:=0 To 2 Do NewXYZ[LoopA]:=XYZ[LoopA];
  End;
End;

// Get Limits From Model
Function GetMappingModelLimits(ModelID: Integer4): Array6Doubles;
Var
  LoopA: Integer4;
  Error: Integer4;
  TotalNodes: Integer4;
  XYZ: Array3Doubles;
Begin
  {Initialize Result}
  {MinX}
  Result[0]:=MaxDouble;
  {MaxX}
  Result[1]:=-MaxDouble;
  {MinY}
  Result[2]:=MaxDouble;
  {MaxY}
  Result[3]:=-MaxDouble;
  {MinZ}
  Result[4]:=MaxDouble;
  {MaxZ}
  Result[5]:=-MaxDouble;

  // Get Total Nodes
  Error:=St7GetTotal(ModelID,tyNODE,TotalNodes);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Get Limits Through Node Searching
  For LoopA:=1 To TotalNodes Do
  Begin
    // Get Node Coords
    Error:=St7GetNodeXYZ(ModelID,LoopA,XYZ);
    If (Error<>ERR7_NoError) Then FatalError(Error);

    // Store Extreme Values
    If (Result[0]>XYZ[0]) Then Result[0]:=XYZ[0];
    If (Result[1]<XYZ[0]) Then Result[1]:=XYZ[0];
    If (Result[2]>XYZ[1]) Then Result[2]:=XYZ[1];
    If (Result[3]<XYZ[1]) Then Result[3]:=XYZ[1];
    If (Result[4]>XYZ[2]) Then Result[4]:=XYZ[2];
    If (Result[5]<XYZ[2]) Then Result[5]:=XYZ[2];
  End;
End;

// Check If Node is Within Limits
Function IsNodeWithinLimits(XYZ: Array3Doubles;MappingLimits: Array6Doubles): Boolean;
Begin
  Result:=(XYZ[0]>MappingLimits[0])And
          (XYZ[0]<MappingLimits[1])And
          (XYZ[1]>MappingLimits[2])And
          (XYZ[1]<MappingLimits[3])And
          (XYZ[2]>MappingLimits[4])And
          (XYZ[2]<MappingLimits[5]);
End;

// Perform Mesh Mapping
Function PerformMeshMapping(MainModelFileName: String;
                            MappingModelFileName: String;
                            MappingModelResultName: String;
                            OutputFileName: String;
                            ParamString: String): Integer4;
Const
  MainModelID = 1;
  MappingModelID = 2;
Var
  LoopA,LoopB,LoopC: Integer4;
  ErrorInt: Integer4;
  TotalStates: Integer4;
  MorphStates: DoubleArray;
  ScratchPathCH: CharString;
  MainModelFileNameCH,MappingModelFileNameCH: CharString;
  MappingModelResultNameCH,SpectralNameCH: CharString;
  OutputFileNameCH: CharString;
  Error: Integer4;
  TotalMainModelNodes: Integer4;
  XYZ: Array3Doubles;
  NewXYZ: Array3Doubles;
  NodeResult: NodeResultArray;
  NumPrimary,NumSecondary: Integer4;
  MappingLimits: Array6Doubles;
Begin
  // Init Result
  Result:=0;

  // Get Parameter Value
  ErrorInt:=ExtractMorphingParam(ParamString,TotalStates,MorphStates);
  If (ErrorInt<>0) Then
  Begin
    Result:=ErrorInt;
    Exit;
  End;

  // Open Main Model
  StrPcopy(MainModelFileNameCH,MainModelFileName);
  StrPcopy(ScratchPathCH,'G:\Temp');
  Error:=St7OpenFile(MainModelID,MainModelFileNameCH,ScratchPathCH);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Open Mapping Model
  StrPcopy(MappingModelFileNameCH,MappingModelFileName);
  StrPcopy(ScratchPathCH,'G:\Temp');
  Error:=St7OpenFile(MappingModelID,MappingModelFileNameCH,ScratchPathCH);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Open Mapping Model Result File
  StrPcopy(MappingModelResultNameCH,MappingModelResultName);
  StrPcopy(SpectralNameCH,'');
  Error:=St7OpenResultFile(MappingModelID,MappingModelResultNameCH,SpectralNameCH,FALSE,NumPrimary,NumSecondary);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Get Mapping Model Limits
  MappingLimits:=GetMappingModelLimits(MappingModelID);

  // Get Main Model Nodes
  Error:=St7GetTotal(MainModelID,tyNODE,TotalMainModelNodes);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Create New Custom Result File
  StrPCopy(OutputFileNameCH,OutputFileName);
  Error:=St7NewResFile(MainModelID,OutputFileNameCH,stLinearBucklingSolver);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Set Custom Result File Paramteres
  Error:=St7SetResFileNumCases(MainModelID,TotalStates);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Set Multiplier for Every Case
  For LoopA:=1 To TotalStates Do
  Begin
    Error:=St7SetResFileMode(MainModelID,LoopA,MorphStates[LoopA]);
    If (Error<>ERR7_NoError) Then FatalError(Error);
  End;

  // Loop Through States
  For LoopA:=1 To TotalStates Do
  Begin
    // Collect Current Parameter
    For LoopB:=1 To TotalMainModelNodes Do
    Begin
      // Create Progress Box
      CreateProgressBox('State ['+IntToStr(LoopA)+'/'+IntToStr(TotalStates)+']; Mapping Node ['+IntToStr(LoopB)+'/'+IntToStr(TotalMainModelNodes)+']...');
      // Get Current Node Coordinates
      Error:=St7GetNodeXYZ(MainModelID,LoopB,XYZ);
      If (Error<>ERR7_NoError) Then FatalError(Error);

      // Check if Nodes if Within Limits
      If (IsNodeWithinLimits(XYZ,MappingLimits)) Then
      Begin
        // Transform Quantity based on Trans Model
        Error:=TransformNode(MappingModelID,XYZ,MorphStates[LoopA],NewXYZ);
        If (Error<>ERR7_NoError) Then FatalError(Error);
      End Else Begin
        // No Diplacement
        For LoopC:=0 To 2 Do NewXYZ[LoopC]:=XYZ[LoopC];
      End;

      // Copy Node Results
      For LoopC:=0 To 5 Do NodeResult[LoopC]:=0.0;
      For LoopC:=0 To 2 Do NodeResult[LoopC]:=NewXYZ[LoopC]-XYZ[LoopC];
      // Set Node Result
      Error:=St7SetResFileNodeResult(MainModelID,LoopA,LoopB,kNodeDisp,NodeResult);
      If (Error<>ERR7_NoError) Then FatalError(Error);
    End;
  End;

  // Get Rid Of Progress Box
  DestroyProgressBox;

  // Close Custom Result File
  Error:=St7CloseResFile(MainModelID);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Close Models
  // Main Model
  Error:=St7CloseFile(MainModelID);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Close Mapping Model Result File
  Error:=St7CloseResultFile(MappingModelID);
  If (Error<>ERR7_NoError) Then FatalError(Error);

  // Trans Model
  Error:=St7CloseFile(MappingModelID);
  If (Error<>ERR7_NoError) Then FatalError(Error);
End;

end.*/

