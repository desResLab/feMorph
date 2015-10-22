unit Sch_Numeric;

interface

Uses SysUtils,General_Var,Dialogs,Math,fft,ap,ssolve;

Const
  MathZero = 1e-8;
  {Errors}
  ip_MCG_NoError = 0;
  ip_MCG_ErrorInPreconditioner = 1;
  ip_MCG_ErrorInSparseMultiply = 2;

Type
  TSchIntegralType = (itSchBase,itSchX);
  TPreconditionerType = (ptJacobi,ptIncompleteCholesky);
  TProbabilityDistributionTypes = (pdtGaussian);
  TRelationTypes = (rtDirect,rtInverse);
  TQuadratureSchemeTypes = (qsClenshawCurtis,qsTrapezoidal);
  TIntegrationPointSupportTypes = (isLegendre{[-1,1]},isHaar{[0,1]});

{Extended Math Functions}
Function  Sign_Sch(Value: Double): Double;
Function  Cube(Value: Double): Double;
Function  EvalFactorial(n: Integer4): Integer4;
{Internal and External Product}
Function  Do3DExternalProduct(V1,V2: Array3Double): Array3Double;
Function  Do3DInternalProduct(V1,V2: Array3Double): Double;OverLoad;
Function  Do3DInternalProduct(V1,V2: Array3Double;NormalizeFirst: Boolean): Double;OverLoad;
Function  Do2DInternalProduct(V1,V2: Array2Double): Double; Overload;
{Normalization}
{Vector}
Procedure Normalize2DVector(Var Vector: Array2Double);
Procedure Normalize3DVector(Var Vector: Array3Double);OverLoad;
Procedure Normalize3DVector(Var Vector: Array3Single);OverLoad;
{Is Null}
Function IsIntegerVectorNull(Size: Integer4;IntVector: Integer4Array): Boolean;
{Matrices}
Procedure TransformToUnitaryColums(RowCount,ColCount: Integer4;OrigAMat: Double2DArray;Var ColNorms: DoubleArray;Var Mat: Double2DArray);
{Scaling}
Procedure Scale3DVector(Var Versor: Array3Double;ScaleFactor: Double);
{Swapping}
Procedure Swap2Values(Var FirstNode,SecondNode: Integer4);
Procedure SwapI4ArrayValues(Index1,Index2: Integer4; Var I4Array: Integer4Array);
Procedure SwapR8ArrayValues(Index1,Index2: Integer4; Var R8Array: DoubleArray);
{Norm Evaluation}
Function  DoEucNorm_0Based(Dim: Integer4;V: DoubleArray): Double;
Function  DoEucNorm_1Based(Dim: Integer4;V: DoubleArray): Double;
Function  DoEucNorm(V: Array3Double): Double;OverLoad;
Function  DoEucNorm(V: Array2Double): Double;OverLoad;
Function  DoEucNorm(V: Array3Single): Double;OverLoad;
Function  EucNorm(Size: Integer4;Vector: DoubleArray): Double;
{2D and 3D Rotations}
{Point in 2D}
Procedure PerformPlaneRotation(Var CurrentXCoord,CurrentYCoord: Double;ShiftX,ShiftY: Double;RotAngle: Double);
{Vectors in 2D and 3D}
Procedure Rotate3DVectorAroundAxis(Var Vector: Array3Double;Angle: Double{In Degrees};Axis: Array3Double);
Function  Rotate2DVector(Vector: Array2Double;Angle: Double{In Degrees}): Array2Double;
{Matrices}
Procedure Rotate2DTensor(Var LocalInertia: Array2x2Double;RotAngle: Double{In Gradi});Overload;
Procedure Rotate2DTensor(LocalInertia: Array2x2Double;RotAngle: Double{In Gradi};Var GlobalInertia: Array2x2Double);Overload;
Procedure TrasformStressTensor(InitialTensor,RotationMatrix: Array3x3Double;Var NewTensor: Array3x3Double);
Function  Trasform3DVector(AMat: Array3x3Double;Vector: Array3Double;Transpose: Boolean): Array3Double;
Procedure Form3DTransformationMatrix(Var AMat: Array3x3Double;V1,V2,V3: Array3Double);
Procedure ExtractSystemAxis(AxisSystem: Array3x3Double;Var Axis1,Axis2,Axis3: Array3Double);
Function  Invert3x3Matrix(Mat: Double2DArray): Double2DArray;
Function  InvertMatrix(Size: Integer4;Mat: Double2DArray): Double2DArray;
Procedure CopyMatrix(Size: Integer4;Source_Mat: Double2DArray;Var Target_Mat: Double2DArray);
Procedure CopyVector(Size: Integer4;Source: DoubleArray;Var Target: DoubleArray);
{Printing Utilities For Matrices}
Procedure PrintMatrixToCSV(MatFileName: String;TotalRows,TotalCols: Integer4;AMat: Double2DArray);
Procedure PrintMatrixToTXT_0Based(MatFileName: String;TotalRows,TotalCols: Integer4;AMat: Double2DArray);
Procedure PrintVectorToTXT_0Based(VecFileName: String;Size: Integer4;Vector: DoubleArray);
{Matrix Algebra}
Function  IsSymmetric(Size: Integer4;AMat: Double2dArray;Var Error: Double): Boolean;
{Full Matrix Multiplication}
Function FullMultiply_0Based(RowCount,ColCount: Integer4;Mat: Double2DArray;Vector: DoubleArray;Var ResultantVector: DoubleArray): String;
Function FullMultiply_T_0Based(RowCount,ColCount: Integer4;Mat: Double2DArray;Vector: DoubleArray;Var ResultantVector: DoubleArray): String;
Function FullMultiply_1Based(RowCount,ColCount: Integer4;Mat: Double2DArray;Vector: DoubleArray;Var ResultantVector: DoubleArray): String;
Function FullMultiply_T_1Based(RowCount,ColCount: Integer4;Mat: Double2DArray;Vector: DoubleArray;Var ResultantVector: DoubleArray): String;
{Sparse Matrix}
Function SparseMultiply(MatSize: Integer4;Diag,Row: Integer4Array;Coeff: DoubleArray; Vector: DoubleArray; Var ResultantVector: DoubleArray): Boolean;
Function AssembleSparseMatrixCoeff(Coeff: Double;Row,Column: Integer4;MatSize: Integer4;MatDiag,MatRow: Integer4Array;var MatCoeff: DoubleArray): Boolean;
Function IterativeMCGSolve(MatSize: Integer4;Mat_Diag,Mat_Row: Integer4Array;Mat_Coeff: DoubleArray;Var RHS: DoubleArray;Preconditioner: TPreconditionerType): Integer4;
{Check Relative Difference Between two Vectors}
Function CheckSolutionDifferece(Size: Integer4;FirstVec: DoubleArray;SecondVec: DoubleArray): Double;
{Linear Equation Sets}
Function GaussElimination(Size: Integer4;AMat: Double2dArray;Var Rhs: DoubleArray;Var ResidualNorm: Double): Integer4;
{Numerical Quadrature}
Function  Eval1DIntegral(Size: Integer4;XValues,YValues: DoubleArray;IntegralType: TSchIntegralType): Double;
{2D and 3D Geometry}
Function  Eval2DPointToSegmentDistance(Px,Py: Double;S1x,S1y,S2x,S2y: Double): Double;Overload;
Function  Eval2DPointToSegmentDistance(Px,Py: Double;S1x,S1y,S2x,S2y: Double; var Param: Double): Double;Overload;
Function  Eval3DPointToSegmentDistance(X0,X1,X2: Array3Double): Double;Overload;
Function  Eval3DPointToSegmentDistance(X0,X1,X2: Array3Double;Var Param: Double): Double;Overload;
Function  Eval3DPointToSegmentProjection(X0,X1,X2: Array3Double): Array3Double;
Function  Eval3DPointToFaceDistance(Px,Py,Pz: Double;T1x,T1y,T1z,T2x,T2y,T2z,T3x,T3y,T3z: Double;Var Distance: Double;Var IsInside: Boolean): Boolean;
Function  Eval2DPointToSegmentVector(Px,Py: Double;S1x,S1y,S2x,S2y: Double): Array2Double;
Function  Eval2DPointVectorSegmentIntersection(N1x,N1y,N2x,N2y,N3x,N3y,N4x,N4y: Double;Var S_Param,T_Param,Ix,Iy: Double): Boolean;
Function  Eval3DPointVectorTriangleIntersection(Ox,Oy,Oz,Vx,Vy,Vz,T1x,T1y,T1z,T2x,T2y,T2z,T3x,T3y,T3z: Double;Var T_Param: Double;Var IsInside: Boolean;Var Ix,Iy,Iz: Double): Boolean;
Function  PointIsInternalToTriangle3D(Point,Vertex1,Vertex2,Vertex3: Array3Double;Var AreaCoord1,AreaCoord2,AreaCoord3: Double): Boolean;
Function  EvalXYTriangleArea(P1,P2,P3: Array3Double): Double;
Procedure Project3DVectorToPlane(Var Vector: Array3Double; PlaneCoeff: Array4Double);
Function  EvalSphereBeamIntersection(SphereCentre: Array3Double;SphereRadius: Double; XYZ1,XYZ2: Array3Double;Var IntersectionParam: Double): Boolean;
{Non linear Equations}
Function SolveIVDegreeEquation(IV_Coeff,III_Coeff,II_Coeff,I_Coeff,Noto_Coeff: Double;Var Final_Value: Double): Boolean;
Function SolveIIIDegreeEquation(III_Coeff,II_Coeff,I_Coeff,Noto_Coeff: Double;Var Final_Value: Double): Boolean;
{Polinomial Evaluation}
Function EvalCubicPolinomal(CurrentCoord: Double;Coeff: Array4Double): Double;
{Table Interpolation}
Function InterpolateTable(RowValue,ColValue: Double;RowCount,ColCount: Integer4;Table: Double2DArray): Double;
{Curve Fitting}
Function PiecewiseInterp(Value: Double;Dim: Integer4;XValues,YValues: DoubleArray): Double;
Function CalculateLeastSquareFit(Dim: Integer4;Degree: Integer;XValues,YValues: DoubleArray;NewX: Double): Double;
Function Eval3DLeastSquaresPlane(Count: Integer4;XValues,YValues,ZValues: DoubleArray;Var PlaneCoeffs: Array4Double): Boolean;
{Probability Distribution}
Function EvalProbabilityCDF(Mean,StDev: Double;Var XValue: Double;CDFType: TProbabilityDistributionTypes;RelationType: TRelationTypes): Double;
Function EvalGaussian(Mean,StDev,XValue: Double): Double;
{Quadrature Formulae: Single}
Function EvalIntegrationPoints(QuadratureRule: TQuadratureSchemeTypes;
                               IntSupport: TIntegrationPointSupportTypes;
                               PointNumber: Integer4;
                               Var CCPoints, CCWeights: DoubleArray): Boolean;
{Quadrature Formulae: Doubled}
Function EvalDoubledIntegrationPoints(QuadratureRule: TQuadratureSchemeTypes;
                                      IntSupport: TIntegrationPointSupportTypes;
                                      Var PointNumber: Integer4;
                                      Var CCPoints, CCWeights: DoubleArray): Boolean;
{Check Small Coefficients}
Procedure CheckSmallCoeffs(ColCount: Integer4;
                           TotalInsertedCols: Integer4;
                           InsertedCols: Integer4Array;
                           ItSol: DoubleArray;
                           ThresholdValue: Double;
                           Var Count: Integer4);

Function EvalClenshawCurtisIntegrationPoints(PointNumber: Integer4;Var CCPoints, CCWeights: DoubleArray): Boolean;

implementation

Uses Sch_TempLists;

Procedure PrintMatrixToTXT_0Based(MatFileName: String;TotalRows,TotalCols: Integer4;AMat: Double2DArray);
Var
  outFile: TextFile;
  LoopA: Integer4;
  LoopB: Integer4;
Begin
  AssignFile(outFile,MatFileName);
  Rewrite(outFile);
  For LoopA:=0 To TotalRows-1 Do
  Begin
    For LoopB:=0 To TotalCols-1 Do
    Begin
      Write(outFile,FloatToStr(AMat[LoopA,LoopB])+' ');
    End;
    Writeln(outFile);
  End;
  CloseFile(outFile);
End;

Procedure PrintVectorToTXT_0Based(VecFileName: String;Size: Integer4;Vector: DoubleArray);
Var
  outFile: TextFile;
  LoopA: Integer4;
  LoopB: Integer4;
Begin
  AssignFile(outFile,VecFileName);
  Rewrite(outFile);
  For LoopA:=0 To Size-1 Do
  Begin
    Writeln(outFile,FloatToStr(Vector[LoopA])+' ');
  End;
  CloseFile(outFile);
End;


Procedure CheckSmallCoeffs(ColCount: Integer4;
                           TotalInsertedCols: Integer4;
                           InsertedCols: Integer4Array;
                           ItSol: DoubleArray;
                           ThresholdValue: Double;
                           Var Count: Integer4);
Var
  LoopA: Integer4;
  Average: Double;
  CurrentCoeff: Double;
Begin
  {Check If There Are Low Coefficients}
  Average:=0.0;
  Count:=0;
  {Eval The Average Of The Coefficients in the Set}
  For LoopA:=1 To ColCount Do
  Begin
    If IsInTempList(LoopA,TotalInsertedCols,InsertedCols) Then
    Begin
      CurrentCoeff:=ABS(ItSol[LoopA]);
      Average:=Average+CurrentCoeff;
    End;
  End;
  Average:=(Average/TotalInsertedCols);
  {Find Small Coeffs}
  For LoopA:=1 To ColCount Do
  Begin
    If IsInTempList(LoopA,TotalInsertedCols,InsertedCols) Then
    Begin
      CurrentCoeff:=ABS(ItSol[LoopA]);
      If (CurrentCoeff<ThresholdValue*Average) Then
      Begin
        Inc(Count);
      End;
    End;
  End;
End;

Function Cube(Value: Double): Double;
Begin
  Result:=Value*Value*Value;
End;

Function EvalProbabilityCDF(Mean,StDev: Double;
                            Var XValue: Double;
                            CDFType: TProbabilityDistributionTypes;
                            RelationType: TRelationTypes): Double;
Begin
  (*Case CDFType Of
    pdtGaussian: Begin
                   If (RelationType=rtDirect) Then
                   Begin
                     Result:=
                   End Else Begin
                     Result:=
                   End;
                 End
  Else
    MessageDlg('Distibution Not Implemented!',mtError,[mbOK],0);
  End;*)
End;

Procedure ProjectOnProbabilityCDF(TotalValues: Integer4;
                                  Var Values: DoubleArray;
                                  CDFType: TProbabilityDistributionTypes;
                                  RelationType: TRelationTypes;
                                  {Distribution Parameters}
                                  Mean,StDev: Double);
Var
  LoopA: Integer4;
Begin
  For LoopA:=1 To TotalValues Do
  Begin
    EvalProbabilityCDF(Mean,StDev,Values[LoopA],CDFType,RelationType);
  End;
End;

function EvalFactorial(n: Integer4): Integer4;
var
  f: LongInt;
  i: Integer;
begin
  f := 1;
  for i := 2 to n do f := f * i;
  Result := f;
end;

Function EvalGaussian(Mean,StDev,XValue: Double): Double;
Begin
  Result:=(1.0/(Sqrt(2.0*PI)*StDev))*(Exp(-((Sqr(XValue-Mean))/(2.0*Sqr(StDev)))));
End;

Function EvalXYTriangleArea(P1,P2,P3: Array3Double): Double;
Var
  First,Second,Third: Double;
Begin
  First:=P2[1]*P3[2]-P3[1]*P2[2];
  Second:=P3[1]*P1[2]-P1[1]*P3[2];
  Third:=P1[1]*P2[2]-P2[1]*P1[2];
  Result:=First+Second+Third;
End;

Function EvalSegno(P0,P1,O: Array3Double): Double;
Begin
  Result:=O[1]*(P0[2]-P1[2])-O[2]*(P0[1]-P1[1])+P0[1]*P1[2]-P0[2]*P1[1];
End;

Function PointIsInternalToTriangle3D(Point,Vertex1,Vertex2,Vertex3: Array3Double;Var AreaCoord1,AreaCoord2,AreaCoord3: Double): Boolean;
Var
  LoopA: Integer4;
  FirstVector,SecondVector,TriangleNormal: Array3Double;
  RotMat: Array3x3Double;
  Segno1,Segno2,Segno3: Double;
  Area1,Area2,Area3,Area: Double;
Begin
  Result:=FALSE;
  {Costruisci la trasformazione matriciale nel piano del Triangolo}
  For LoopA:=1 To 3 Do FirstVector[LoopA]:=Vertex2[LoopA]-Vertex1[LoopA];
  Normalize3DVector(FirstVector);
  For LoopA:=1 To 3 Do SecondVector[LoopA]:=Vertex3[LoopA]-Vertex2[LoopA];
  Normalize3DVector(SecondVector);
  TriangleNormal:=Do3DExternalProduct(FirstVector,SecondVector);
  Normalize3DVector(TriangleNormal);
  SecondVector:=Do3DExternalProduct(TriangleNormal,FirstVector);
  Normalize3DVector(SecondVector);
  For LoopA:=1 To 3 Do RotMat[LoopA,1]:=FirstVector[LoopA];
  For LoopA:=1 To 3 Do RotMat[LoopA,2]:=SecondVector[LoopA];
  For LoopA:=1 To 3 Do RotMat[LoopA,3]:=TriangleNormal[LoopA];
  {Trasforma Tutti i Punti nel Piano}
  Point:=Trasform3DVector(RotMat,Point,TRUE);
  Vertex1:=Trasform3DVector(RotMat,Vertex1,TRUE);
  Vertex2:=Trasform3DVector(RotMat,Vertex2,TRUE);
  Vertex3:=Trasform3DVector(RotMat,Vertex3,TRUE);
  {Valuta Le Funzioni Segno}
  Segno1:=EvalSegno(Vertex1,Vertex2,Vertex3)*EvalSegno(Vertex1,Vertex2,Point);
  Segno2:=EvalSegno(Vertex2,Vertex3,Vertex1)*EvalSegno(Vertex2,Vertex3,Point);
  Segno3:=EvalSegno(Vertex3,Vertex1,Vertex2)*EvalSegno(Vertex3,Vertex1,Point);
  If (Segno1<-MathZero)Or(Segno2<-MathZero)Or(Segno3<-MathZero) Then
  Begin
    Result:=FALSE;
    AreaCoord1:=0.0;
    AreaCoord2:=0.0;
    AreaCoord3:=0.0;
  End Else Begin
    Result:=TRUE;
    {Calcolo Delle Coordinate d'area}
    Area1:=EvalXYTriangleArea(Point,Vertex2,Vertex3);
    Area2:=EvalXYTriangleArea(Point,Vertex3,Vertex1);
    Area3:=EvalXYTriangleArea(Point,Vertex1,Vertex2);
    Area:=EvalXYTriangleArea(Vertex1,Vertex2,Vertex3);
    {Calcolo Coord}
    AreaCoord1:=(Area1/Area);
    AreaCoord2:=(Area2/Area);
    AreaCoord3:=(Area3/Area);
  End;
End;

Function EvalPreconditioner(MatSize: Integer4;
                            Mat_Diag,Mat_Row: Integer4Array;
                            Mat_Coeff: DoubleArray;
                            Var P_Diag,P_Row: Integer4Array;
                            Var P_Coeff: DoubleArray;
                            Preconditioner: TPreconditionerType): Boolean;
Var
  LoopA: Integer4;
Begin
  Result:=TRUE;
  {Valuta il precondizionatore}
  Case Preconditioner Of
    ptJacobi: Begin
                {Alloca le Matrici del Precondizionatore}
                SetLength(P_Diag,MatSize+1);
                SetLength(P_Row,MatSize+1);
                SetLength(P_Coeff,MatSize+1);
                {Riempi le Matrici}
                For LoopA:=1 To MatSize Do
                Begin
                  P_Diag[LoopA]:=LoopA;
                  P_Row[LoopA]:=LoopA;
                  If ABS(Mat_Coeff[Mat_Diag[LoopA]])<MathZero Then
                  Begin
                    Result:=FALSE;
                    Exit;
                  End;
                  P_Coeff[LoopA]:=1.0/(Mat_Coeff[Mat_Diag[LoopA]]);
                End;
              End;
    ptIncompleteCholesky: Begin
                            Result:=FALSE;
                          End;
  End;
End;

{Risoluzione con il gradiente coniugato precondizionato}
Function IterativeMCGSolve(MatSize: Integer4;Mat_Diag,Mat_Row: Integer4Array;Mat_Coeff: DoubleArray;Var RHS: DoubleArray;Preconditioner: TPreconditionerType): Integer4;
Const
  IterationTolerance = 1e-12;
  IterationMax = 10000;
Var
  LoopA: Integer4;
  Scarto: Double;
  CurrentIteration: Integer4;
  P_Diag,P_Row: Integer4Array;
  P_Coeff: DoubleArray;
  {Vettori Iterazioni}
  X_k,R_k,P_k,AP_k: DoubleArray;
  X_kp1,R_kp1,P_kp1,PR_kp1,PAP_k: DoubleArray;
  {Quantità Scalari}
  Alfa_k_Num,Alfa_k_Den: Double;
  Beta_k_Num: Double;
  Alfa_k,Beta_k: Double;
Begin
  Result:=ip_MCG_NoError;
  {Allocazione dei vettori di Iterazione}
  SetLength( X_k,MatSize+1);
  SetLength(R_k,MatSize+1);
  SetLength(P_k,MatSize+1);
  SetLength(AP_k,MatSize+1);
  SetLength(X_kp1,MatSize+1);
  SetLength(R_kp1,MatSize+1);
  SetLength(P_kp1,MatSize+1);
  SetLength(PR_kp1,MatSize+1);
  SetLength(PAP_k,MatSize+1);
  {Calcola il Precondizionatore}
  If Not(EvalPreconditioner(MatSize,Mat_Diag,Mat_Row,Mat_Coeff,P_Diag,P_Row,P_Coeff,ptJacobi)) Then
  Begin
    Result:=ip_MCG_ErrorInPreconditioner;
  End;
  {Determina la Prima Iterata}
  If Not(SparseMultiply(MatSize,P_Diag,P_Row,P_Coeff,RHS,X_k)) Then
  Begin
    Result:=ip_MCG_ErrorInSparseMultiply;
  End;
  {Determina Il Residuo}
  If Not(SparseMultiply(MatSize,Mat_Diag,Mat_Row,Mat_Coeff,X_k,R_k)) Then
  Begin
    Result:=ip_MCG_ErrorInSparseMultiply;
  End;
  For LoopA:=1 To MatSize Do R_k[LoopA]:=RHS[LoopA]-R_k[LoopA];
  {Calcola lo scarto}
  Scarto:=0.0;
  For LoopA:=1 To MatSize Do Scarto:=Scarto+R_k[LoopA]*R_k[LoopA];
  Scarto:=Sqrt(Scarto);
  {Calcola P0 o pk}
  If Not(SparseMultiply(MatSize,P_Diag,P_Row,P_Coeff,R_k,P_k)) Then
  Begin
    Result:=ip_MCG_ErrorInSparseMultiply;
  End;
  CurrentIteration:=0;
  {Schema Iterativo}
  While (Scarto>IterationTolerance)And(CurrentIteration<IterationMax) Do
  Begin
    Inc(CurrentIteration);
    {Calcola Apk}
    If Not(SparseMultiply(MatSize,Mat_Diag,Mat_Row,Mat_Coeff,p_k,AP_k)) Then
    Begin
      Result:=ip_MCG_ErrorInSparseMultiply;
    End;
    Alfa_k_Num:=0.0;
    Alfa_k_Den:=0.0;
    For LoopA:=1 To MatSize Do
    Begin
      Alfa_k_Num:=Alfa_k_Num+P_k[LoopA]*R_k[LoopA];
      Alfa_k_Den:=Alfa_k_Den+P_k[LoopA]*AP_k[LoopA];
    End;
    {Calcola Alfak}
    Alfa_k:=(Alfa_k_Num/Alfa_k_Den);
    {Calcola X_kp1}
    For LoopA:=1 To MatSize Do
    Begin
      X_kp1[LoopA]:=X_k[LoopA]+Alfa_k*P_k[LoopA];
    End;
    {Calcola Rkp1}
    For LoopA:=1 To MatSize Do
    Begin
      R_kp1[LoopA]:=R_k[LoopA]-Alfa_k*AP_k[LoopA];
    End;
    {Calcola Prkp1}
    If Not(SparseMultiply(MatSize,P_Diag,P_Row,P_Coeff,R_kp1,PR_kp1)) Then
    Begin
      Result:=ip_MCG_ErrorInSparseMultiply;
    End;
    {Calcola PAPk}
    If Not(SparseMultiply(MatSize,P_Diag,P_Row,P_Coeff,AP_k,PAP_k)) Then
    Begin
      Result:=ip_MCG_ErrorInSparseMultiply;
    End;
    {Calcola Betak}
    Beta_k_Num:=0.0;
    For LoopA:=1 To MatSize Do
    Begin
      Beta_k_Num:=Beta_k_Num+R_kp1[LoopA]*PAP_k[LoopA];
    End;
    Beta_k:=-(Beta_k_Num/Alfa_k_Den);
    {Calcola Pkp1}
    For LoopA:=1 To MatSize Do
    Begin
      P_kp1[LoopA]:=PR_kp1[LoopA]+Beta_k*P_k[LoopA];
    End;
    {Calcola Scarto}
    Scarto:=0.0;
    For LoopA:=1 To MatSize Do Scarto:=Scarto+R_kp1[LoopA]*R_kp1[LoopA];
    Scarto:=Sqrt(Scarto);
    {Aggiorna Variabili di Ciclo}
    For LoopA:=1 To MatSize Do
    Begin
      X_k[LoopA]:=X_kp1[LoopA];
      P_k[LoopA]:=P_kp1[LoopA];
      R_k[LoopA]:=R_kp1[LoopA];
    End;
  End;
  {Ricopia la soluzione nel RHS}
  For LoopA:=1 To MatSize Do RHS[LoopA]:=X_k[LoopA];
  {Deallocate Preconditioner}
  SetLength(P_Diag,0);
  SetLength(P_Row,0);
  SetLength(P_Coeff,0);
  FreeMemory(P_Diag);
  FreeMemory(P_Row);
  FreeMemory(P_Coeff);
  {Deallocazione dei vettori di Iterazione}
  SetLength(X_k,0);
  SetLength(R_k,0);
  SetLength(P_k,0);
  SetLength(AP_k,0);
  SetLength(X_kp1,0);
  SetLength(R_kp1,0);
  SetLength(P_kp1,0);
  SetLength(PR_kp1,0);
  SetLength(PAP_k,0);
  FreeMemory(X_k);
  FreeMemory(R_k);
  FreeMemory(P_k);
  FreeMemory(AP_k);
  FreeMemory(X_kp1);
  FreeMemory(R_kp1);
  FreeMemory(P_kp1);
  FreeMemory(PR_kp1);
  FreeMemory(PAP_k);
End;

{La matrice di rigidezza memorizza le colonne della Upper diagonal (Ovviamente matrice simmetrica)}
Function SparseMultiply(MatSize: Integer4;Diag,Row: Integer4Array;Coeff: DoubleArray; Vector: DoubleArray; Var ResultantVector: DoubleArray): Boolean;
Var
  LoopA,LoopB: Integer4;
  StartIndex,FinishIndex: Integer4;
Begin
  Result:=TRUE;
  {Initializza la componente del vettore risultante}
  For LoopA:=1 To MatSize Do
  Begin
    ResultantVector[LoopA]:=0.0;
  End;
  {Moltiplicazione}
  For LoopA:=1 To MatSize Do
  Begin
    {Scegli l'indice di partenza}
    If LoopA=1 Then StartIndex:=1
    Else StartIndex:=(Diag[LoopA-1]+1);
    FinishIndex:=Diag[LoopA];
    {Prima Parte}
    For LoopB:=StartIndex To FinishIndex Do
    Begin
      {Prima Parte}
      ResultantVector[LoopA]:=ResultantVector[LoopA]+Coeff[LoopB]*Vector[Row[LoopB]];
      {Seconda Parte}
      If Row[LoopB]<>LoopA Then ResultantVector[Row[LoopB]]:=ResultantVector[Row[LoopB]]+Coeff[LoopB]*Vector[LoopA];      
    End;
  End;
End;

Function EvalCubicPolinomal(CurrentCoord: Double;Coeff: Array4Double): Double;
Var
  LoopA: Integer4;
Begin
  Result:=0.0;
  For LoopA:=1 To 4 Do
  Begin
    Result:=Result+Power(CurrentCoord,4-LoopA)*Coeff[LoopA];
  End;
End;

Function IsSymmetric(Size: Integer4;AMat: Double2dArray;Var Error: Double): Boolean;
Var
  LoopA,LoopB: Integer4;
Begin
  Error:=0.0;
  Result:=TRUE;
  For LoopA:=1 To Size Do
  Begin
    For LoopB:=LoopA+1 To Size Do
    Begin
      If (ABS(AMat[LoopA,LoopB]-AMat[LoopB,LoopA])>Error) Then Error:=ABS(AMat[LoopA,LoopB]-AMat[LoopB,LoopA]);
      If (ABS(AMat[LoopA,LoopB]-AMat[LoopB,LoopA])>MathZero) Then
      Begin
        Result:=FALSE;
        Exit;
      End;
    End;
  End;
End;

Function Eval1DIntegral(Size: Integer4;XValues,YValues: DoubleArray;IntegralType: TSchIntegralType): Double;
Var
  LoopA: Integer4;
  SommaBasi: Double;
Begin
  Result:=0.0;
  SommaBasi:=0.0;
  For LoopA:=1 To Size-1 Do
  Begin
    Case IntegralType Of
      itSchBase: SommaBasi:=(YValues[LoopA+1]+YValues[LoopA]);
      itSchX: SommaBasi:=(YValues[LoopA+1]*XValues[LoopA+1]+YValues[LoopA]*XValues[LoopA])
    End;
    Result:=Result+0.5*SommaBasi*(XValues[LoopA+1]-XValues[LoopA]);
  End;
End; 

Function Do3DExternalProduct(V1,V2: Array3Double): Array3Double;
Begin
  Result[1]:=V1[2]*V2[3]-V2[2]*V1[3];
  Result[2]:=V1[3]*V2[1]-V2[3]*V1[1];
  Result[3]:=V1[1]*V2[2]-V2[1]*V1[2];
End;

Function DoEucNorm(V: Array3Double): Double;OverLoad;
Var
  A: Integer4;
  Sum: Double;
Begin
  Sum:=0.0;
  For A:=1 To 3 Do
  Begin
    Sum:=Sum+Sqr(V[A]);
  End;
  Result:=Sqrt(Sum);
End;

Function DoEucNorm(V: Array3Single): Double;OverLoad;
Var
  A: Integer4;
  Sum: Double;
Begin
  Sum:=0.0;
  For A:=1 To 3 Do
  Begin
    Sum:=Sum+Sqr(V[A]);
  End;
  Result:=Sqrt(Sum);
End;


Function DoEucNorm(V: Array2Double): Double;OverLoad;
Var
  A: Integer4;
  Sum: Double;
Begin
  Sum:=0.0;
  For A:=1 To 2 Do
  Begin
    Sum:=Sum+Sqr(V[A]);
  End;
  Result:=Sqrt(Sum);
End;


Function DoEucNorm_1Based(Dim: Integer4;V: DoubleArray): Double;OverLoad;
Var
  A: Integer4;
  Sum: Double;
Begin
  Sum:=0.0;
  For A:=1 To Dim Do
  Begin
    Sum:=Sum+Sqr(V[A]);
  End;
  Result:=Sqrt(Sum);
End;

Function DoEucNorm_0Based(Dim: Integer4;V: DoubleArray): Double;OverLoad;
Var
  A: Integer4;
  Sum: Double;
Begin
  Sum:=0.0;
  For A:=0 To Dim-1 Do
  Begin
    Sum:=Sum+Sqr(V[A]);
  End;
  Result:=Sqrt(Sum);
End;


Procedure SortCurveArray(Dim: Integer4;Var XValues,YValues: DoubleArray);
Var
  LoopA,LoopB: Integer4;
  FirstValue,SecondValue: Double;
  Dum: Double;
Begin
  For LoopA:=1 To Dim Do
  Begin
    For LoopB:=(LoopA+1) To Dim Do
    Begin
      FirstValue:=XValues[LoopA];
      SecondValue:=XValues[LoopB];
      If (FirstValue>SecondValue) Then
      Begin
        {Swap X And Y}
        {X}
        Dum:=XValues[LoopA];
        XValues[LoopA]:=XValues[LoopB];
        XValues[LoopB]:=Dum;
        {Y}
        Dum:=YValues[LoopA];
        YValues[LoopA]:=YValues[LoopB];
        YValues[LoopB]:=Dum;
      End;
    End;
  End;
End;

Function PiecewiseInterp(Value: Double;Dim: Integer4;XValues,YValues: DoubleArray): Double;
Var
  Count: Integer4;
  Found: Boolean;
  StartX,StartY,FinishX,FinishY: Double;

Begin
  Count:=1;
  Found:=FALSE;
  Result:=0.0;
  {Sort Vectors for Increasing X}
  SortCurveArray(Dim,XValues,YValues);
  While (Count<Dim)And(Not(Found)) Do
  Begin
    Inc(Count);
    StartX:=XValues[Count-1];
    StartY:=YValues[Count-1];
    FinishX:=XValues[Count];
    FinishY:=YValues[Count];
    If (Value>=StartX)And(Value<=FinishX) Then
    Begin
      Found:=TRUE;
      Result:=StartY+((FinishY-StartY)/(FinishX-StartX))*(Value-StartX);
    End;
  End;
  If Not(Found) Then
  Begin
    Try
    If (Value>XValues[Dim]) Then Result:=YValues[Dim]
    Else If (Value<XValues[1]) Then Result:=YValues[1];
    Except
      ShowMessage('T');
    End;

  End;
End;

Function EucNorm(Size: Integer4;Vector: DoubleArray): Double;
Var
  A: Integer4;
Begin
  Result:=0.0;
  For A:=1 To Size Do
  Begin
    Result:=Result+Sqr(Vector[A]);
  End;
  Result:=Sqrt(Result);
End;

Function GaussElimination(Size: Integer4;AMat: Double2dArray;Var Rhs: DoubleArray;Var ResidualNorm: Double): Integer4;
Var
  A,B,C: Integer4;
  InitialNorm: Double;
  Factor: Double;
  Solution,Residual: DoubleArray;
Begin
  Result:=0;
  Factor:=0.0;
  {Allocate}
  SetLength(Solution,Size+1);
  SetLength(Residual,Size+1);
  InitialNorm:=EucNorm(Size,Rhs);
  // Reduction To Triangular Form
  For A:=1 To Size Do
  Begin
    For B:=(A+1) To Size Do
    Begin
      If ABS(AMat[A,A])<MathZero Then
      Begin
        Result:=1;
        Exit;
      End Else Begin
        Factor:=-(AMat[B,A]/AMat[A,A]);
      End;
      For C:=1 To Size Do
      Begin
        AMat[B,C]:=AMat[B,C]+Factor*AMat[A,C];
      End;
      // Rhs
      Rhs[B]:=Rhs[B]+Factor*Rhs[A];
    End;
  End;
  // Back Sub
  For A:=Size DownTo 1 Do
  Begin
    Solution[A]:=Rhs[A];
    For B:=A+1 To Size Do
    Begin
      Solution[A]:=Solution[A]-AMat[A,B]*Solution[B];
    End;
    If ABS(AMat[A,A])<MathZero Then
      Begin
        Result:=1;
        Exit;
      End Else Begin
        Solution[A]:=Solution[A]/AMat[A,A];
      End;
  End;
  // Verify Error
  For A:=1 To Size Do
  Begin
    Residual[A]:=Rhs[A];
    For B:=1 To Size Do
    Begin
      Residual[A]:=Residual[A]-AMat[A,B]*Solution[B];
    End;
  End;
  // Copy the Solution Vector
  For A:=1 To Size Do
  Begin
    Rhs[A]:=Solution[A];
  End;
  If (ABS(InitialNorm)<MathZero) Then ResidualNorm:=0.0
  Else ResidualNorm:=(EucNorm(Size,Residual)/InitialNorm);
  {De-Allocate}
  SetLength(Solution,0);
  SetLength(Residual,0);
  FreeMemory(Solution);
  FreeMemory(Residual);
End;

{Retta approssimante ai minimi quadrati}
Function CalculateLeastSquareFit(Dim: Integer4;Degree: Integer;XValues,YValues: DoubleArray;NewX: Double): Double;
Var
  A,B,C: Integer4;
  s1,v0,s2,v1,t: Double;
  MediaX,MediaY: Double;
  A1,B1,Sum: Double;
  AMat: Double2dArray;
  QMat: Double2dArray;
  RHS: DoubleArray;
  ResidualNorm: Double;
Begin
  If Degree=1 Then
  Begin
    {Calcolo delle medie}
    s1:=0.0;
    v0:=0.0;
    s2:=0.0;
    v1:=0.0;
    t:=0.0;
    For A:=1 To Dim Do
    Begin
      s1:=s1+XValues[A];
      v0:=v0+YValues[A];
      s2:=s2+Sqr(XValues[A]);
      v1:=v1+XValues[A]*YValues[A];
      t:=t+Sqr(YValues[A]);
    End;
    MediaX:=(s1/Dim);
    MediaY:=(v0/Dim);
    A1:=s2-Dim*Sqr(MediaX);
    B1:=v1-Dim*MediaX*MediaY;
    Result:=(B1/A1)*NewX+(MediaY-(B1/A1)*MediaX);
  End Else Begin
    {Calcolo della matrice A}
    {Allocazione}
    Setlength(AMat,Dim+1);
    For A:=1 To Dim Do Setlength(AMat[A],Degree+2);
    For A:=1 To Dim Do
    Begin
      For B:=1 To Degree+1 Do
      Begin
        AMat[A,B]:=Power(XValues[A],B-1);
      End;
    End;
    {Determinazione della matrice QMat}
    SetLength(QMat,Degree+2);
    For A:=1 To Degree+1 Do SetLength(QMat[A],Degree+2);
    For A:=1 To Degree+1 Do
    Begin
      For B:=1 To Degree+1 Do
      Begin
        Sum:=0.0;
        For C:=1 To Dim Do
        Begin
          Sum:=Sum+AMat[C,A]*AMat[C,B];
        End;
        QMat[A,B]:=Sum;
      End;
    End;
    {Deallocate AMat}
    For A:=1 To Dim Do
    Begin
      Setlength(AMat[A],0);
      FreeMemory(AMat[A]);
    End;
    Setlength(AMat,0);
    FreeMemory(AMat);
    {Calculate the RHS}
    {Allocate}
    SetLength(RHS,Degree+2);
    For A:=1 To Degree+1 Do
    Begin
      Sum:=0.0;
      For B:=1 To Dim Do
      Begin
        Sum:=Sum+YValues[B]*Power(XValues[B],A-1);
      End;
      RHS[A]:=Sum;
    End;
    {Solve the Linear Set of equations}
    If (GaussElimination(Degree+1,QMat,RHS,ResidualNorm)<>0) Then
    Begin
      MessageDlg('Errore nell''inversione della matrice ai minimi quadrati!',mtError,[mbOk],0);
    End;
    {Write Result}
    Result:=0.0;
    For A:=1 To Degree+1 Do
    Begin
      Result:=Result+RHS[A]*Power(NewX,A-1);
    End;
    {Deallocate RHS e QMat}
    For A:=1 To Degree+1 Do
    Begin
      Setlength(QMat[A],0);
      FreeMemory(QMat[A]);
    End;
    Setlength(QMat,0);
    FreeMemory(QMat);
  End;
End;

Procedure Normalize2DVector(Var Vector: Array2Double);
Var
  Modulus: Double;
Begin
  Modulus:=Sqrt(Sqr(Vector[1])+Sqr(Vector[2]));
  If (Modulus>0.0) Then Vector[1]:=Vector[1]/Modulus;
  If (Modulus>0.0) Then Vector[2]:=Vector[2]/Modulus;
End;

Procedure Normalize3DVector(Var Vector: Array3Double);Overload;
Var
  Modulus: Double;
Begin
  Modulus:=Sqrt(Sqr(Vector[1])+Sqr(Vector[2])+Sqr(Vector[3]));
  If (Modulus>0.0) Then Vector[1]:=Vector[1]/Modulus;
  If (Modulus>0.0) Then Vector[2]:=Vector[2]/Modulus;
  If (Modulus>0.0) Then Vector[3]:=Vector[3]/Modulus;
End;

Procedure Normalize3DVector(Var Vector: Array3Single);Overload;
Var
  Modulus: Double;
Begin
  Modulus:=Sqrt(Sqr(Vector[1])+Sqr(Vector[2])+Sqr(Vector[3]));
  If (Modulus>0.0) Then Vector[1]:=Vector[1]/Modulus;
  If (Modulus>0.0) Then Vector[2]:=Vector[2]/Modulus;
  If (Modulus>0.0) Then Vector[3]:=Vector[3]/Modulus;
End;


Procedure Rotate3DVectorAroundAxis(Var Vector: Array3Double;Angle: Double{Angolo in Gradi};Axis: Array3Double);
Var
  LoopA,LoopB: Integer4;
  w,x,y,z,Sum: Double;
  QuatMat: Array4x4Double;
  QuatVect,QuatRes: Array4Double;
Begin
  {Normalizza il vettore}
  Normalize3DVector(Axis);
  {Scrivi le componenti del quaternione}
  w:=cos(0.5*Angle*(Pi/180.0));
  x:=Axis[1]*sin(0.5*Angle*(Pi/180.0));
  y:=Axis[2]*sin(0.5*Angle*(Pi/180.0));
  z:=Axis[3]*sin(0.5*Angle*(Pi/180.0));
  {Trova la matrice di Rotazione}
  QuatMat[1,1]:=Sqr(w)+Sqr(x)-Sqr(y)-Sqr(z);
  QuatMat[1,2]:=2.0*x*y+2.0*w*z;
  QuatMat[1,3]:=2.0*x*z-2.0*w*y;
  QuatMat[1,4]:=0.0;
  QuatMat[2,1]:=2.0*x*y-2.0*w*z;
  QuatMat[2,2]:=Sqr(w)-Sqr(x)+Sqr(y)-Sqr(z);
  QuatMat[2,3]:=2.0*y*z+2.0*w*x;
  QuatMat[2,4]:=0.0;
  QuatMat[3,1]:=2.0*x*z+2.0*w*y;
  QuatMat[3,2]:=2.0*y*z-2.0*w*x;
  QuatMat[3,3]:=Sqr(w)-Sqr(x)-Sqr(y)+Sqr(z);
  QuatMat[3,4]:=0.0;
  QuatMat[4,1]:=0.0;
  QuatMat[4,2]:=0.0;
  QuatMat[4,3]:=0.0;
  QuatMat[4,4]:=Sqr(w)+Sqr(x)+Sqr(y)+Sqr(z);
  {Espandi il vettore da ruotare}
  QuatVect[4]:=0.0;
  For LoopA:=1 To 3 Do QuatVect[LoopA]:=Vector[LoopA];
  {Esegui la moltiplicazione matriciale}
  For LoopA:=1 To 4 Do
  Begin
    Sum:=0.0;
    For LoopB:=1 To 4 Do
    Begin
      Sum:=Sum+QuatVect[LoopB]*QuatMat[LoopB,LoopA];
    End;
    QuatRes[LoopA]:=Sum;
  End;
  {Copia il risultato finale}
  For LoopA:=1 To 3 Do Vector[LoopA]:=QuatRes[LoopA];
End;

Function  Rotate2DVector(Vector: Array2Double;Angle: Double{In Degrees}): Array2Double;
Var
  RadAngle: Double;
Begin
  RadAngle:=Angle*(Pi/180.0);
  Result[1]:=Vector[1]*Cos(RadAngle)+Vector[2]*Sin(RadAngle);
  Result[2]:=-Vector[1]*Sin(RadAngle)+Vector[2]*Cos(RadAngle);
End;

Function SolveIVDegreeEquation(IV_Coeff,III_Coeff,II_Coeff,I_Coeff,Noto_Coeff: Double;Var Final_Value: Double): Boolean;
Function EvalCurrentFunction(XValue: Double): Double;
Begin
  Result:=IV_Coeff*Power(XValue,4.0)+III_Coeff*Power(XValue,3.0)+II_Coeff*Sqr(XValue)+I_Coeff*XValue+Noto_Coeff;
End;
Function EvalCurrentDerivative(XValue: Double): Double;
Begin
  Result:=4.0*IV_Coeff*Power(XValue,3.0)+3.0*III_Coeff*Sqr(XValue)+2.0*II_Coeff*XValue+I_Coeff;
End;
Const
  IterationTolerance = 1e-5;
  MaxIterations = 1000;
Var
  CurrentXValue,CurrentNorm: Double;
  IterationNumber: Integer4;
  CurrentFunction,CurrentDerivative: Double;
  NewXValue: Double;
Begin
  Result:=TRUE;
  CurrentXValue:=0.0;
  CurrentNorm:=10.0;
  IterationNumber:=1;
  While (CurrentNorm>IterationTolerance)And(IterationNumber<MaxIterations) Do
  Begin
    CurrentFunction:=EvalCurrentFunction(CurrentXValue);
    CurrentDerivative:=EvalCurrentDerivative(CurrentXValue);
    NewXValue:=CurrentXValue-(CurrentFunction/CurrentDerivative);
    If (IterationNumber>1) Then CurrentNorm:=Abs((NewXValue-CurrentXValue)/CurrentXValue)
    Else CurrentNorm:=10.0;
    CurrentXValue:=NewXValue;
    {Update}
    Inc(IterationNumber);
  End;
  If IterationNumber>MaxIterations Then
  Begin
    Result:=FALSE;
    Exit;
  End;
  Final_Value:=CurrentXValue;
End;

Function SolveIIIDegreeEquation(III_Coeff,II_Coeff,I_Coeff,Noto_Coeff: Double;Var Final_Value: Double): Boolean;
Function EvalCurrentFunction(XValue: Double): Double;
Begin
  Result:=III_Coeff*Power(XValue,3.0)+II_Coeff*Sqr(XValue)+I_Coeff*XValue+Noto_Coeff;
End;
Function EvalCurrentDerivative(XValue: Double): Double;
Begin
  Result:=3.0*III_Coeff*Sqr(XValue)+2.0*II_Coeff*XValue+I_Coeff;
End;
Const
  IterationTolerance = 1e-5;
  MaxIterations = 1000;
Var
  CurrentXValue,CurrentNorm: Double;
  IterationNumber: Integer4;
  CurrentFunction,CurrentDerivative: Double;
  NewXValue: Double;
Begin
  Result:=TRUE;
  CurrentXValue:=0.0;
  CurrentNorm:=10.0;
  IterationNumber:=1;
  While (CurrentNorm>IterationTolerance)And(IterationNumber<MaxIterations) Do
  Begin
    CurrentFunction:=EvalCurrentFunction(CurrentXValue);
    CurrentDerivative:=EvalCurrentDerivative(CurrentXValue);
    NewXValue:=CurrentXValue-(CurrentFunction/CurrentDerivative);
    If (IterationNumber>1) Then CurrentNorm:=Abs((NewXValue-CurrentXValue)/CurrentXValue)
    Else CurrentNorm:=10.0;
    CurrentXValue:=NewXValue;
    {Update}
    Inc(IterationNumber);
  End;
  If IterationNumber>MaxIterations Then
  Begin
    Result:=FALSE;
    Exit;
  End;
  Final_Value:=CurrentXValue;
End;

Procedure PerformPlaneRotation(Var CurrentXCoord,CurrentYCoord: Double;ShiftX,ShiftY: Double;RotAngle: Double{In Gradi});
Var
  NewXCoord,NewYCoord: Double;
Begin
  NewXCoord:=(CurrentXCoord-ShiftX)*Cos(RotAngle*(Pi/180.0))-(CurrentYCoord-ShiftY)*Sin(RotAngle*(Pi/180.0));
  NewYCoord:=(CurrentXCoord-ShiftX)*Sin(RotAngle*(Pi/180.0))+(CurrentYCoord-ShiftY)*Cos(RotAngle*(Pi/180.0));
  CurrentXCoord:=NewXCoord;
  CurrentYCoord:=NewYCoord;
End;

Procedure Rotate2DTensor(LocalInertia: Array2x2Double;RotAngle: Double{In Gradi};Var GlobalInertia: Array2x2Double);Overload;
Var
  LoopA,LoopB,LoopC: Integer4;
  RotMat,FirstMat: Array2x2Double;
  Sum: Double;
Begin
  {Determina la matrice di rotazione piana}
  RotMat[1,1]:=Cos(RotAngle*(Pi/180.0));
  RotMat[1,2]:=-Sin(RotAngle*(Pi/180.0));
  RotMat[2,1]:=Sin(RotAngle*(Pi/180.0));
  RotMat[2,2]:=Cos(RotAngle*(Pi/180.0));
  {Esegui la trasformazione}
  For LoopA:=1 To 2 Do
  Begin
    For LoopB:=1 To 2 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 2 Do
      Begin
        Sum:=Sum+RotMat[LoopC,LoopA]*LocalInertia[LoopC,LoopB];
      End;
      FirstMat[LoopA,LoopB]:=Sum;
    End;
  End;
  For LoopA:=1 To 2 Do
  Begin
    For LoopB:=1 To 2 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 2 Do
      Begin
        Sum:=Sum+FirstMat[LoopA,LoopC]*RotMat[LoopC,LoopB];
      End;
      GlobalInertia[LoopA,LoopB]:=Sum;
    End;
  End;
End;

Procedure Rotate2DTensor(Var LocalInertia: Array2x2Double;RotAngle: Double{In Gradi});
Var
  LoopA,LoopB,LoopC: Integer4;
  RotMat,FirstMat: Array2x2Double;
  Sum: Double;
Begin
  {Determina la matrice di rotazione piana}
  RotMat[1,1]:=Cos(RotAngle*(Pi/180.0));
  RotMat[1,2]:=-Sin(RotAngle*(Pi/180.0));
  RotMat[2,1]:=Sin(RotAngle*(Pi/180.0));
  RotMat[2,2]:=Cos(RotAngle*(Pi/180.0));
  {Esegui la trasformazione}
  For LoopA:=1 To 2 Do
  Begin
    For LoopB:=1 To 2 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 2 Do
      Begin
        Sum:=Sum+RotMat[LoopC,LoopA]*LocalInertia[LoopC,LoopB];
      End;
      FirstMat[LoopA,LoopB]:=Sum;
    End;
  End;
  For LoopA:=1 To 2 Do
  Begin
    For LoopB:=1 To 2 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 2 Do
      Begin
        Sum:=Sum+FirstMat[LoopA,LoopC]*RotMat[LoopC,LoopB];
      End;
      LocalInertia[LoopA,LoopB]:=Sum;
    End;
  End;
End;

Function  Do3DInternalProduct(V1,V2: Array3Double): Double; OverLoad;
Var
  LoopA: Integer2;
Begin
  Result:=0.0;
  For LoopA:=1 To 3 Do Result:=Result+V1[LoopA]*V2[LoopA];
End;

Function Do3DInternalProduct(V1,V2: Array3Double;NormalizeFirst: Boolean): Double; Overload;
Var
  LoopA: Integer2;
Begin
  IF NormalizeFirst Then
  Begin
    Normalize3DVector(V1);
    Normalize3DVector(V2);
  End;
  Result:=0.0;
  For LoopA:=1 To 3 Do Result:=Result+V1[LoopA]*V2[LoopA];
End;

Function Do2DInternalProduct(V1,V2: Array2Double): Double; Overload;
Var
  LoopA: Integer2;
Begin
  Result:=0.0;
  For LoopA:=1 To 2 Do Result:=Result+V1[LoopA]*V2[LoopA];
End;

Procedure TrasformStressTensor(InitialTensor,RotationMatrix: Array3x3Double;Var NewTensor: Array3x3Double);
Var
  LoopA,LoopB,LoopC: Integer4;
  Sum: Double;
  TempMat: Array3x3Double;
Begin
  {Pre-moltiplicazione}
  For LoopA:=1 To 3 Do
  Begin
    For LoopB:=1 To 3 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 3 Do
      Begin
        Sum:=Sum+RotationMatrix[LoopC,LoopA]*InitialTensor[LoopC,LoopB]
      End;
      TempMat[LoopA,LoopB]:=Sum;
    End;
  End;
  {Post-moltiplicazione}
  For LoopA:=1 To 3 Do
  Begin
    For LoopB:=1 To 3 Do
    Begin
      Sum:=0.0;
      For LoopC:=1 To 3 Do
      Begin
        Sum:=Sum+TempMat[LoopA,LoopC]*RotationMatrix[LoopC,LoopB]
      End;
      NewTensor[LoopA,LoopB]:=Sum;
    End;
  End;
End;

Function Trasform3DVector(AMat: Array3x3Double;Vector: Array3Double;Transpose: Boolean): Array3Double;
Var
  LoopA,LoopB: Integer2;
Begin
  For LoopA:=1 To 3 Do
  Begin
    Result[LoopA]:=0.0;
    For LoopB:=1 To 3 Do
    Begin
      If Transpose Then Result[LoopA]:=Result[LoopA]+AMat[LoopB,LoopA]*Vector[LoopB]
      Else Result[LoopA]:=Result[LoopA]+AMat[LoopA,LoopB]*Vector[LoopB];
    End;
  End;
End;

Procedure Form3DTransformationMatrix(Var AMat: Array3x3Double;V1,V2,V3: Array3Double);
Var
  LoopA: Integer2;
Begin
  Normalize3DVector(V1);
  Normalize3DVector(V2);
  Normalize3DVector(V3);
  For LoopA:=1 To 3 Do AMat[LoopA,1]:=V1[LoopA];
  For LoopA:=1 To 3 Do AMat[LoopA,2]:=V2[LoopA];
  For LoopA:=1 To 3 Do AMat[LoopA,3]:=V3[LoopA];
End;

Function Sign_Sch(Value: Double): Double;
Begin
  Result:=Sign(Value);
  If ABS(Result)<MathZero Then Result:=1.0;
End;

Function Eval3DLeastSquaresPlane(Count: Integer4;XValues,YValues,ZValues: DoubleArray;Var PlaneCoeffs: Array4Double): Boolean;
Var
  LoopA,LoopB: Integer4;
  AMat: Double2DArray;
  Rhs: DoubleArray;
  ResidualNorm: Double;
Begin
  Result:=TRUE;
  {Alloca ed Initializa la matrice dei coefficienti}
  SetLength(AMat,3+1);
  SetLength(Rhs,3+1);
  For LoopA:=1 To 3 Do
  Begin
    SetLength(AMat[LoopA],3+1);
  End;
  {Initializza}
  For LoopA:=1 To 3 Do
  Begin
    Rhs[LoopA]:=0.0;
    For LoopB:=1 To 3 Do
    Begin
      AMat[LoopA,LoopB]:=0.0;
    End;
  End;
  {Determina i coefficienti della matrice}
  For LoopA:=2 To Count Do
  Begin
    {Prima Riga}
    AMat[1,1]:=AMat[1,1]+(XValues[Count]-XValues[1])*(XValues[Count]-XValues[1]);
    AMat[1,2]:=AMat[1,2]+(XValues[Count]-XValues[1])*(YValues[Count]-YValues[1]);
    AMat[1,3]:=AMat[1,3]+(XValues[Count]-XValues[1])*(ZValues[Count]-ZValues[1]);
    {Seconda Riga}
    AMat[2,1]:=AMat[2,1]+(YValues[Count]-YValues[1])*(XValues[Count]-XValues[1]);
    AMat[2,2]:=AMat[2,2]+(YValues[Count]-YValues[1])*(YValues[Count]-YValues[1]);
    AMat[2,3]:=AMat[2,3]+(YValues[Count]-YValues[1])*(ZValues[Count]-ZValues[1]);
    {Terza Riga}
    AMat[3,1]:=AMat[3,1]+(ZValues[Count]-ZValues[1])*(XValues[Count]-XValues[1]);
    AMat[3,2]:=AMat[3,2]+(ZValues[Count]-ZValues[1])*(YValues[Count]-YValues[1]);
    AMat[3,3]:=AMat[3,3]+(ZValues[Count]-ZValues[1])*(ZValues[Count]-ZValues[1]);
    {Noto}
    Rhs[1]:=Rhs[1]-(XValues[Count]-XValues[1]);
    Rhs[2]:=Rhs[2]-(YValues[Count]-YValues[1]);
    Rhs[3]:=Rhs[3]-(ZValues[Count]-ZValues[1]);
  End;
  {Solve Equation Set}
  If (GaussElimination(3,AMat,Rhs,ResidualNorm)=0) Then
  Begin
    For LoopA:=1 To 3 Do PlaneCoeffs[LoopA]:=Rhs[LoopA];
    PlaneCoeffs[4]:=1.0;
  End Else Begin
    Result:=FALSE;
    Exit;
  End;
  {Deallocate}
  For LoopA:=1 To 3 Do
  Begin
    SetLength(AMat[LoopA],0);
  End;
  SetLength(AMat,0);
  SetLength(Rhs,0);
  FreeMemory(AMat);
  FreeMemory(Rhs);
End;

Procedure Project3DVectorToPlane(Var Vector: Array3Double; PlaneCoeff: Array4Double);
Var
  LoopA: Integer4;
  NormalVersor: Array3Double;
  NormalComponent: Double;
Begin
  {Eval Plane Versor}
  For LoopA:=1 To 3 Do NormalVersor[LoopA]:=PlaneCoeff[LoopA];
  Normalize3DVector(NormalVersor);
  {Find the Projection}
  NormalComponent:=Do3DInternalProduct(Vector,NormalVersor);
  For LoopA:=1 To 3 Do Vector[LoopA]:=Vector[LoopA]-NormalComponent*NormalVersor[LoopA];
End;

Function  EvalSphereBeamIntersection(SphereCentre: Array3Double;SphereRadius: Double; XYZ1,XYZ2: Array3Double;Var IntersectionParam: Double): Boolean;
Var
  Error: Integer4;
  FirstNodeIncluded,SecondNodeIncluded: Boolean;
  A_Coeff,B_Coeff,C_Coeff,Discriminant: Double;
  SegmentParam,SegmentParam1: Double;
Begin
  {Eval State}
  FirstNodeIncluded:=(Sqr(XYZ1[1]-SphereCentre[1])+Sqr(XYZ1[2]-SphereCentre[2])+Sqr(XYZ1[3]-SphereCentre[3])<Sqr(SphereRadius));
  SecondNodeIncluded:=(Sqr(XYZ2[1]-SphereCentre[1])+Sqr(XYZ2[2]-SphereCentre[2])+Sqr(XYZ2[3]-SphereCentre[3])<Sqr(SphereRadius));
  Result:=((FirstNodeIncluded)Xor(SecondNodeIncluded));
  If Not(Result) Then
  Begin
    IntersectionParam:=0.0;
    Exit;
  End;
  {Eval Coeffs}
  A_Coeff:=Sqr(XYZ2[1]-XYZ1[1])+Sqr(XYZ2[2]-XYZ1[2])+Sqr(XYZ2[3]-XYZ1[3]);
  B_Coeff:=2.0*((XYZ2[1]-XYZ1[1])*(XYZ1[1]-SphereCentre[1])+(XYZ2[2]-XYZ1[2])*(XYZ1[2]-SphereCentre[2])+(XYZ2[3]-XYZ1[3])*(XYZ1[3]-SphereCentre[3]));
  C_Coeff:=Sqr(SphereCentre[1])+Sqr(SphereCentre[2])+Sqr(SphereCentre[3])+Sqr(XYZ1[1])+Sqr(XYZ1[2])+Sqr(XYZ1[3])- 2.0*(SphereCentre[1]*XYZ1[1]+SphereCentre[2]*XYZ1[2]+SphereCentre[3]*XYZ1[3])-Sqr(SphereRadius);
  Discriminant:=Sqr(B_Coeff)-4.0*A_Coeff*C_Coeff;
  If (Discriminant<-MathZero) Then
  Begin
    ShowMessage('Errore Interno: Intersezione Beam Sfera');
  End Else Begin
    SegmentParam:=(-B_Coeff+Sqrt(Discriminant))/(2.0*A_Coeff);
    SegmentParam1:=(-B_Coeff-Sqrt(Discriminant))/(2.0*A_Coeff);
    If (SegmentParam>0.0-MathZero)And(SegmentParam<1.0+MathZero) Then
    Begin
      IntersectionParam:=SegmentParam;
    End Else If (SegmentParam1>0.0-MathZero)And(SegmentParam1<1.0+MathZero) Then
    Begin
      IntersectionParam:=SegmentParam1;
    End Else Begin
      ShowMessage('Errore Interno: Parametro Errato: '+FloatToStrF(SegmentParam,ffFixed,15,5));
    End;
  End;
End;

Procedure ExtractSystemAxis(AxisSystem: Array3x3Double;Var Axis1,Axis2,Axis3: Array3Double);
Var
  LoopA: Integer2;
Begin
  {Recupera le Colonne}
  For LoopA:=1 To 3 Do
  Begin
    Axis1[LoopA]:=AxisSystem[LoopA,1];
    Axis2[LoopA]:=AxisSystem[LoopA,2];
    Axis3[LoopA]:=AxisSystem[LoopA,3];
  End;
End;

{Assemble a coefficient in a Sparse Matrix}
Function AssembleSparseMatrixCoeff(Coeff: Double;Row,Column: Integer4;MatSize: Integer4;MatDiag,MatRow: Integer4Array;var MatCoeff: DoubleArray): Boolean;
Var
  LoopA: Integer4;
  StartIndex,EndIndex: Integer4;
Begin
  Result:=FALSE;
  If (Row=1) Then StartIndex:=1
  Else StartIndex:=MatDiag[Row-1]+1;
  EndIndex:=(MatDiag[Row]);
  For LoopA:=StartIndex To EndIndex Do
  Begin
    If (MatRow[LoopA]=Column) Then
    Begin
      MatCoeff[LoopA]:=MatCoeff[LoopA]+Coeff;
      Result:=TRUE;
    End;
  End;
End;

Function  Eval2DPointToSegmentDistance(Px,Py: Double;S1x,S1y,S2x,S2y: Double): Double;
Var
  SegDist,UValue: Double;
  Intx,Inty: Double;
Begin
   SegDist:=Sqr(S2x-S1x)+Sqr(S2y-S1y);
   UValue:=(((Px-S1x)*(S2x-S1x)+(Py-S1y)*(S2y-S1y))/SegDist);
   Intx:=S1x+(S2x-S1x)*UValue;
   Inty:=S1y+(S2y-S1y)*UValue;
   Result:=Sqrt(Sqr(Px-Intx)+Sqr(Py-Inty));
End;

Function Eval2DPointToSegmentDistance(Px,Py: Double;S1x,S1y,S2x,S2y: Double; var Param: Double): Double;Overload;
Var
  SegDist: Double;
  Intx,Inty: Double;
Begin
   SegDist:=Sqr(S2x-S1x)+Sqr(S2y-S1y);
   Param:=(((Px-S1x)*(S2x-S1x)+(Py-S1y)*(S2y-S1y))/SegDist);
   Intx:=S1x+(S2x-S1x)*Param;
   Inty:=S1y+(S2y-S1y)*Param;
   Result:=Sqrt(Sqr(Px-Intx)+Sqr(Py-Inty));
End;

Function  Eval3DPointToSegmentDistance(X0,X1,X2: Array3Double): Double;Overload;
Var
  LoopA: Integer4;
  Diff1,Diff2,IntPoint: Array3Double;
  Numerator,Denominator,TParam: Double;
Begin
   For LoopA:=1 To 3 Do Diff1[LoopA]:=X1[LoopA]-X0[LoopA];
   For LoopA:=1 To 3 Do Diff2[LoopA]:=X2[LoopA]-X1[LoopA];

   Numerator:=Do3DInternalProduct(Diff1,Diff2);
   Denominator:=Sqr(DoEucNorm(Diff2));

   TParam:=-(Numerator)/(Denominator);

   For LoopA:=1 To 3 Do IntPoint[LoopA]:=X1[LoopA]+(X2[LoopA]-X1[LoopA])*TParam;

   Result:=Sqrt(Sqr(IntPoint[1]-X0[1])+Sqr(IntPoint[2]-X0[2])+Sqr(IntPoint[3]-X0[3]));
End;

Function  Eval3DPointToSegmentDistance(X0,X1,X2: Array3Double;Var Param: Double): Double;Overload;
Var
  LoopA: Integer4;
  Diff1,Diff2,IntPoint: Array3Double;
  Numerator,Denominator: Double;
Begin
   For LoopA:=1 To 3 Do Diff1[LoopA]:=X1[LoopA]-X0[LoopA];
   For LoopA:=1 To 3 Do Diff2[LoopA]:=X2[LoopA]-X1[LoopA];

   Numerator:=Do3DInternalProduct(Diff1,Diff2);
   Denominator:=Sqr(DoEucNorm(Diff2));

   Param:=-(Numerator)/(Denominator);

   For LoopA:=1 To 3 Do IntPoint[LoopA]:=X1[LoopA]+(X2[LoopA]-X1[LoopA])*Param;

   Result:=Sqrt(Sqr(IntPoint[1]-X0[1])+Sqr(IntPoint[2]-X0[2])+Sqr(IntPoint[3]-X0[3]));
End;


Function  Eval3DPointToSegmentProjection(X0,X1,X2: Array3Double): Array3Double;
Var
  LoopA: Integer4;
  Diff1,Diff2,IntPoint: Array3Double;
  Numerator,Denominator,TParam: Double;
Begin
   For LoopA:=1 To 3 Do Diff1[LoopA]:=X1[LoopA]-X0[LoopA];
   For LoopA:=1 To 3 Do Diff2[LoopA]:=X2[LoopA]-X1[LoopA];

   Numerator:=Do3DInternalProduct(Diff1,Diff2);
   Denominator:=Sqr(DoEucNorm(Diff2));

   TParam:=-(Numerator)/(Denominator);

   For LoopA:=1 To 3 Do Result[LoopA]:=X1[LoopA]+(X2[LoopA]-X1[LoopA])*TParam;
End;


Function  Eval2DPointToSegmentVector(Px,Py: Double;S1x,S1y,S2x,S2y: Double): Array2Double;
Var
  SegDist,UValue: Double;
  Intx,Inty: Double;
Begin
   SegDist:=Sqr(S2x-S1x)+Sqr(S2y-S1y);
   UValue:=(((Px-S1x)*(S2x-S1x)+(Py-S1y)*(S2y-S1y))/SegDist);
   Intx:=S1x+(S2x-S1x)*UValue;
   Inty:=S1y+(S2y-S1y)*UValue;
   Result[1]:=Intx-Px;
   Result[2]:=Inty-Py;
End;

Function Eval2DPointVectorSegmentIntersection(N1x,N1y,N2x,N2y,N3x,N3y,N4x,N4y: Double;Var S_Param,T_Param,Ix,Iy: Double): Boolean;
Var
  Alfa1,Alfa2,Alfa3: Double;
  Alfa4,Alfa5,Alfa6: Double;
  T_Denom,T_Num: Double;
Begin
  Alfa1:=N2x-N1x;
  Alfa2:=N4x-N3x;
  Alfa3:=N3x-N1x;
  Alfa4:=N2y-N1y;
  Alfa5:=N4y-N3y;
  Alfa6:=N3y-N1y;
  {Find Intersection}
  If ABS(Alfa1)<MathZero Then
  Begin
    If ABS(Alfa2)>MathZero Then
    Begin
      T_Param:=-(Alfa3/Alfa2);
      S_Param:=(Alfa6+T_Param*Alfa5)/(Alfa4);
      Ix:=N1x+S_Param*(N2x-N1x);
      Iy:=N1y+S_Param*(N2y-N1y);
      If (T_Param>-MathZero)And(T_Param<1.0+MathZero) Then
      Begin
        If (S_Param>-MathZero) Then Result:=TRUE
        Else Result:=FALSE;
      End Else Result:=FALSE;
    End Else Begin
      Result:=FALSE;
      S_Param:=0.0;
      T_Param:=0.0;
    End;
  End Else Begin
    T_Denom:=Alfa5-Alfa2*(Alfa4/Alfa1);
    If ABS(T_Denom)<MathZero Then
    Begin
      Result:=FALSE;
      S_Param:=0.0;
      T_Param:=0.0;
    End Else Begin
      T_Num:=Alfa6-Alfa3*(Alfa4/Alfa1);
      T_Param:=-(T_Num)/(T_Denom);
      S_Param:=(Alfa3+T_Param*Alfa2)/(Alfa1);
      Ix:=N1x+S_Param*(N2x-N1x);
      Iy:=N1y+S_Param*(N2y-N1y);
      If (T_Param>-MathZero)And(T_Param<1.0+MathZero) Then
      Begin
        If (S_Param>-MathZero) Then Result:=TRUE
        Else Result:=FALSE;
      End Else Result:=FALSE;
    End;
  End;
End;

Procedure Swap2Values(Var FirstNode,SecondNode: Integer4);
Var
  Dum: Integer4;
Begin
  Dum:=FirstNode;
  FirstNode:=SecondNode;
  SecondNode:=Dum;
End;

Function Eval3DPointVectorTriangleIntersection(Ox,Oy,Oz,
                                               Vx,Vy,Vz,
                                               T1x,T1y,T1z,
                                               T2x,T2y,T2z,
                                               T3x,T3y,T3z: Double;
                                               Var T_Param: Double;
                                               Var IsInside: Boolean;
                                               Var Ix,Iy,Iz: Double): Boolean;
Var
  Aux1,Aux2: Array3Double;
  PlaneNormal: Array3Double;
  Num,Den: Double;
  FirstCheck,SecondCheck,ThirdCheck: Double;
Begin
  {Find The Normal to the plane}
  {1}
  Aux1[1]:=T2x-T1x;
  Aux1[2]:=T2y-T1y;
  Aux1[3]:=T2z-T1z;
  {2}
  Aux2[1]:=T3x-T1x;
  Aux2[2]:=T3y-T1y;
  Aux2[3]:=T3z-T1z;
  {Cross Product}
  PlaneNormal:=Do3DExternalProduct(Aux1,Aux2);
  Normalize3DVector(PlaneNormal);
  Aux1[1]:=Ox-T1x;
  Aux1[2]:=Oy-T1y;
  Aux1[3]:=Oz-T1z;
  {V}
  Aux2[1]:=Vx;
  Aux2[2]:=Vy;
  Aux2[3]:=Vz;
  {Make Internal Products}
  Num:=-Do3DInternalProduct(Aux1,PlaneNormal);
  Den:=Do3DInternalProduct(Aux2,PlaneNormal);
  If ABS(Den)<MathZero Then
  Begin
    Result:=FALSE;
    T_Param:=0.0;
    IsInside:=FALSE;
    Ix:=0.0;
    Iy:=0.0;
    Iz:=0.0;
    Exit;
  End;
  {T_Param}
  T_Param:=(Num)/(Den);
  {Write Intersection}
  Ix:=Ox+T_Param*Vx;
  Iy:=Oy+T_Param*Vy;
  Iz:=Oz+T_Param*Vz;
  {Eval If Point is Inside}
  {Check 01}
  Aux1[1]:=T2x-T1x;
  Aux1[2]:=T2y-T1y;
  Aux1[3]:=T2z-T1z;
  Aux2[1]:=Ix-T1x;
  Aux2[2]:=Iy-T1y;
  Aux2[3]:=Iz-T1z;
  FirstCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check 02}
  Aux1[1]:=T3x-T2x;
  Aux1[2]:=T3y-T2y;
  Aux1[3]:=T3z-T2z;
  Aux2[1]:=Ix-T2x;
  Aux2[2]:=Iy-T2y;
  Aux2[3]:=Iz-T2z;
  SecondCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check 03}
  Aux1[1]:=T1x-T3x;
  Aux1[2]:=T1y-T3y;
  Aux1[3]:=T1z-T3z;
  Aux2[1]:=Ix-T3x;
  Aux2[2]:=Iy-T3y;
  Aux2[3]:=Iz-T3z;
  ThirdCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check If Inside}
  IsInside:=(FirstCheck>-Mathzero)And(SecondCheck>-Mathzero)And(ThirdCheck>-Mathzero);
  If (IsInside)And(T_Param>-MathZero) Then Result:=TRUE
  Else Result:=FALSE;
End;

Function Eval3DPointToFaceDistance(Px,Py,Pz: Double;T1x,T1y,T1z,T2x,T2y,T2z,T3x,T3y,T3z: Double;Var Distance: Double;Var IsInside: Boolean): Boolean;
Var
  Aux1,Aux2: Array3Double;
  PlaneNormal: Array3Double;
  Ix,Iy,Iz: Double;
  FirstCheck,SecondCheck,ThirdCheck: Double;
Begin
  {Find The Normal to the plane}
  {1}
  Aux1[1]:=T2x-T1x;
  Aux1[2]:=T2y-T1y;
  Aux1[3]:=T2z-T1z;
  {2}
  Aux2[1]:=T3x-T1x;
  Aux2[2]:=T3y-T1y;
  Aux2[3]:=T3z-T1z;
  {Cross Product}
  PlaneNormal:=Do3DExternalProduct(Aux1,Aux2);
  Normalize3DVector(PlaneNormal);
  Aux1[1]:=Px-T1x;
  Aux1[2]:=Py-T1y;
  Aux1[3]:=Pz-T1z;
  {Make Internal Products}
  Distance:=ABS(Do3DInternalProduct(Aux1,PlaneNormal));
  {Write Intersection}
  Ix:=Px+Distance*PlaneNormal[1];
  Iy:=Py+Distance*PlaneNormal[2];
  Iz:=Pz+Distance*PlaneNormal[3];
  {Eval If Point is Inside}
  {Check 01}
  Aux1[1]:=T2x-T1x;
  Aux1[2]:=T2y-T1y;
  Aux1[3]:=T2z-T1z;
  Aux2[1]:=Ix-T1x;
  Aux2[2]:=Iy-T1y;
  Aux2[3]:=Iz-T1z;
  FirstCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check 02}
  Aux1[1]:=T3x-T2x;
  Aux1[2]:=T3y-T2y;
  Aux1[3]:=T3z-T2z;
  Aux2[1]:=Ix-T2x;
  Aux2[2]:=Iy-T2y;
  Aux2[3]:=Iz-T2z;
  SecondCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check 03}
  Aux1[1]:=T1x-T3x;
  Aux1[2]:=T1y-T3y;
  Aux1[3]:=T1z-T3z;
  Aux2[1]:=Ix-T3x;
  Aux2[2]:=Iy-T3y;
  Aux2[3]:=Iz-T3z;
  ThirdCheck:=Do3DInternalProduct(Do3DExternalProduct(Aux1,Aux2),PlaneNormal);
  {Check If Inside}
  IsInside:=(FirstCheck>-Mathzero)And(SecondCheck>-Mathzero)And(ThirdCheck>-Mathzero);
  If (IsInside) Then Result:=TRUE
  Else Result:=FALSE;
End;

Function Invert3x3Matrix(Mat: Double2DArray): Double2DArray;
Var
  LoopA: Integer4;
  DetJ: Double;
Begin
  {Valuta il Determinante}
  DetJ:=Mat[1,1]*(Mat[2,2]*Mat[3,3]-Mat[2,3]*Mat[3,2])-
        Mat[1,2]*(Mat[2,1]*Mat[3,3]-Mat[2,3]*Mat[3,1])+
        Mat[1,3]*(Mat[2,1]*Mat[3,2]-Mat[2,2]*Mat[3,1]);
  {Allocate Result}
  SetLength(Result,4);
  For LoopA:=1 To 3 Do SetLength(Result[LoopA],4);
  {Write Matrix}
  Result[1,1]:= (1.0/DetJ)*(Mat[2,2]*Mat[3,3]-Mat[2,3]*Mat[3,2]);
  Result[2,2]:= (1.0/DetJ)*(Mat[1,1]*Mat[3,3]-Mat[1,3]*Mat[3,1]);
  Result[3,3]:= (1.0/DetJ)*(Mat[1,1]*Mat[2,2]-Mat[1,2]*Mat[2,1]);
  Result[1,2]:=-(1.0/DetJ)*(Mat[1,2]*Mat[3,3]-Mat[1,3]*Mat[3,2]);
  Result[1,3]:= (1.0/DetJ)*(Mat[1,2]*Mat[2,3]-Mat[2,2]*Mat[1,3]);
  Result[2,1]:=-(1.0/DetJ)*(Mat[2,1]*Mat[3,3]-Mat[2,3]*Mat[3,1]);
  Result[2,3]:=-(1.0/DetJ)*(Mat[1,1]*Mat[2,3]-Mat[1,3]*Mat[2,1]);
  Result[3,1]:= (1.0/DetJ)*(Mat[2,1]*Mat[3,2]-Mat[3,1]*Mat[2,2]);
  Result[3,2]:=-(1.0/DetJ)*(Mat[1,1]*Mat[3,2]-Mat[1,2]*Mat[3,1]);
End;

Procedure CopyMatrix(Size: Integer4;Source_Mat: Double2DArray;Var Target_Mat: Double2DArray);
Var
  LoopA,LoopB: Integer4;
Begin
  For LoopA:=1 To Size Do
  Begin
    For LoopB:=1 To Size Do
    Begin
      Target_Mat[LoopA,LoopB]:=Source_Mat[LoopA,LoopB];
    End;
  End;
End;

Procedure CopyVector(Size: Integer4;Source: DoubleArray;Var Target: DoubleArray);
Var
  LoopA: Integer4;
Begin
  For LoopA:=1 To Size Do Target[LoopA]:=Source[LoopA];
End;

Function InvertMatrix(Size: Integer4;Mat: Double2DArray): Double2DArray;
Var
  LoopA,LoopB: Integer4;
  RHS: DoubleArray;
  GaussError: Integer4;
  ResidualNorm: Double;
Begin
  SetLength(Result,Size+1);
  For LoopA:=1 To Size Do SetLength(Result[LoopA],Size+1);
  {Allocate RHS}
  SetLength(RHS,Size+1);
  {Loop On the Matrix Colums}
  For LoopA:=1 To Size Do
  Begin
    {Set The RHS}
    For LoopB:=1 To Size Do
    Begin
      If LoopB=LoopA Then RHS[LoopB]:=1.0
      Else RHS[LoopB]:=0.0;
    End;
    {Solve The System}
    GaussError:=GaussElimination(Size,Mat,RHS,ResidualNorm);
    If (GaussError>0) Then ShowMessage('Internal: Could Not Find Solution in InvertMatrix');
    {Copy Column Into Result}
    For LoopB:=1 To Size Do
    Begin
      Result[LoopB,LoopA]:=RHS[LoopB];
    End;
  End;
  {Deallocate RHS}
  SetLength(RHS,0);
  FreeMemory(RHS);
End;

Procedure SwapI4ArrayValues(Index1,Index2: Integer4; Var I4Array: Integer4Array);
Var
  DumI4: Integer4;
Begin
  DumI4:=I4Array[Index1];
  I4Array[Index1]:=I4Array[Index2];
  I4Array[Index2]:=DumI4;
End;
Procedure SwapR8ArrayValues(Index1,Index2: Integer4; Var R8Array: DoubleArray);
Var
  DumR8: Double;
Begin
  DumR8:=R8Array[Index1];
  R8Array[Index1]:=R8Array[Index2];
  R8Array[Index2]:=DumR8;
End;

{Eval Newton-Cotes Intergrals: Trapezoidal Rule [-1;1]}
Function EvalTrapezoidalIntegrationPoints(PointNumber: Integer4;Var CCPoints, CCWeights: DoubleArray): Boolean;
Var
  LoopA: Integer4;
  IntervalValue: Double;
  CurrentPoint: Double;
Begin
  Result:=TRUE;
  If (PointNumber<1) Then
  Begin
    Result:=FALSE;
    Exit;
  End;
  If (PointNumber=1) Then
  Begin
    SetLength(CCPoints,PointNumber+1);
    SetLength(CCWeights,PointNumber+1);
    CCPoints[PointNumber]:=0.0;
    CCWeights[PointNumber]:=2.0;
    Exit;
  End;
  SetLength(CCPoints,PointNumber+1);
  SetLength(CCWeights,PointNumber+1);
  {Single Interval}
  IntervalValue:=2.0/(PointNumber-1);
  {Build Points}
  CurrentPoint:=-1.0;
  For LoopA:=1 To PointNumber Do
  Begin
    CCPoints[LoopA]:=CurrentPoint;
    CurrentPoint:=CurrentPoint+IntervalValue;
    If (LoopA=1)Or(LoopA=PointNumber) Then CCWeights[LoopA]:=0.5*IntervalValue
    Else CCWeights[LoopA]:=IntervalValue
  End;
End;

{Eval Clenshaw-Curtis Integrals}
Function EvalClenshawCurtisIntegrationPoints(PointNumber: Integer4;Var CCPoints, CCWeights: DoubleArray): Boolean;
const
  printDebug = true;
Var
  LoopA: Integer4;
  LoopEnd: Integer4;
  V: DoubleArray;
  G: TComplex1DArray;
  Weights: TReal1DArray;
  w0: Double;
  N_Value: Integer4;
  vString: String;
Begin
  Result:=TRUE;
  If (PointNumber<1) Then
  Begin
    Result:=FALSE;
    Exit;
  End;
  If (PointNumber=1) Then
  Begin
    SetLength(CCPoints,PointNumber+1);
    SetLength(CCWeights,PointNumber+1);
    CCPoints[PointNumber]:=0.0;
    CCWeights[PointNumber]:=2.0;
    Exit;
  End;
  N_Value:=(PointNumber-1);
  LoopEnd:=Trunc(N_Value/2);
  {Allocation}
  SetLength(CCPoints,PointNumber+1);
  SetLength(CCWeights,PointNumber+1);
  SetLength(V,LoopEnd+1);
  SetLength(G,LoopEnd+1);
  For LoopA:=1 To PointNumber Do
  Begin
    CCPoints[LoopA]:=-cos(((LoopA-1)*Pi)/(PointNumber-1))
  End;
  {Form Vector V}
  For LoopA:=0 To LoopEnd Do
  Begin
    If (LoopA=LoopEnd) Then
    Begin
      V[LoopA]:=((N_Value-3.0)/(2.0*((N_Value)Div(2))-1.0))-1.0;
    End Else Begin
      V[LoopA]:=2.0/(1-4.0*Sqr(LoopA));
    End;
  End;
  w0:=(1.0/(Sqr(N_Value)-1.0+((N_Value)Mod(2))));
  {Form Vector G}
  For LoopA:=0 To LoopEnd Do
  Begin
    If (LoopA=LoopEnd) Then
    Begin
      G[LoopA].X:=V[LoopA]+w0*((2.0-((N_Value)Mod(2)))*N_Value-1.0);
      G[LoopA].Y:=0.0;
    End Else Begin
      G[LoopA].X:=V[LoopA]-w0;
      G[LoopA].Y:=0.0;
    End;
  End;

  Try
    {Eval the Inverse DFT to Find the Weights}
    FFTR1DInv(G,PointNumber-1,Weights);

    {Print G}
    (*vString:='Transform'+#13;
    If printDebug Then
    Begin
      For LoopA:=0 To PointNumber-2 Do
      Begin
        vString:=vString+FloaTtoStr(Weights[LoopA])+#13;
      End;
    End;
    ShowMessage(vString);*)

    For LoopA:=1 To PointNumber-1 Do
    Begin
      CCWeights[LoopA]:=Weights[LoopA-1];
      (*If (LoopA<=((PointNumber)Div(2))+1) Then
      Begin
        CCWeights[LoopA]:=Weights[LoopA-1];
      End Else Begin
        CCWeights[LoopA]:=Weights[PointNumber-LoopA];
      End;*)
    End;
    CCWeights[PointNumber]:=Weights[0];
    {Deallocate}
    SetLength(V,0);
    SetLength(G,0);
    FreeMemory(V);
    FreeMemory(G);
  Except
    Result:=FALSE;
    {Deallocate}
    SetLength(V,0);
    SetLength(G,0);
    FreeMemory(V);
    FreeMemory(G);
    Exit;
  End;
End;

{Scale Vector}
Procedure Scale3DVector(Var Versor: Array3Double;ScaleFactor: Double);
Var
  LoopA: Integer4;
Begin
  For LoopA:=1 To 3 Do Versor[LoopA]:=Versor[LoopA]*ScaleFactor;
End;

{General Integration Point Extraction Framework}
Function EvalIntegrationPoints(QuadratureRule: TQuadratureSchemeTypes;
                               IntSupport: TIntegrationPointSupportTypes;
                               PointNumber: Integer4;
                               Var CCPoints, CCWeights: DoubleArray): Boolean;
Var
  LoopA,LoopB: Integer4;
  CurrentValue: Double;
  {Doubled Quadrature}
  TotIntPoints: Integer4;
  IntPoints,IntWeights: DoubleArray;
Begin
  Result:=TRUE;
  If (PointNumber>0) Then
  Begin
    Case QuadratureRule Of
      qsClenshawCurtis: Begin
                          Result:=EvalClenshawCurtisIntegrationPoints(PointNumber,CCPoints,CCWeights);
                        End;
      qsTrapezoidal: Begin
                       Result:=EvalTrapezoidalIntegrationPoints(PointNumber,CCPoints,CCWeights);
                     End;
    End;
    Case IntSupport Of
      isHaar: Begin
                {Scale Points to [0,1]}
                For LoopA:=1 To PointNumber Do
                Begin
                  CurrentValue:=CCPoints[LoopA];
                  CurrentValue:=0.5*(CurrentValue+1);
                  CCPoints[LoopA]:=CurrentValue;
                  CCWeights[LoopA]:=0.5*CCWeights[LoopA];
                End;
              End;
    End;
  End;
End;

{Eval Doubled Integration Rule}
Function EvalDoubledIntegrationPoints(QuadratureRule: TQuadratureSchemeTypes;
                                      IntSupport: TIntegrationPointSupportTypes;
                                      Var PointNumber: Integer4;
                                      Var CCPoints, CCWeights: DoubleArray): Boolean;
Var
  LoopA,LoopB: Integer4;
  CurrentValue: Double;
  {Doubled Quadrature}
  TotIntPoints: Integer4;
  IntPoints,IntWeights: DoubleArray;
Begin
  Result:=TRUE;
  Case QuadratureRule Of
    qsClenshawCurtis: Begin
                        Result:=EvalClenshawCurtisIntegrationPoints(PointNumber,CCPoints,CCWeights);
                      End;
    qsTrapezoidal: Begin
                     Result:=EvalTrapezoidalIntegrationPoints(PointNumber,CCPoints,CCWeights);
                   End;
  End;
  Case IntSupport Of
    isHaar: Begin
              {Scale Points to [0,1]}
              For LoopA:=1 To PointNumber Do
              Begin
                CurrentValue:=CCPoints[LoopA];
                CurrentValue:=0.5*(CurrentValue+1);
                CCPoints[LoopA]:=CurrentValue;
                CCWeights[LoopA]:=0.5*CCWeights[LoopA];
              End;
            End;
  End;

  {Build Double Quadrature Rule}
  TotIntPoints:=2*PointNumber;
  SetLength(IntPoints,TotIntPoints+1);
  SetLength(IntWeights,TotIntPoints+1);

  {Transform Points And Weights}
  For LoopA:=1 To 2 Do
  Begin
    For LoopB:=1 To PointNumber Do
    Begin
      {Map Points}
      If (IntSupport=isHaar) Then IntPoints[(LoopA-1)*PointNumber+LoopB]:=(LoopA-1)*0.5+0.5*CCPoints[LoopB]
      Else IntPoints[(LoopA-1)*PointNumber+LoopB]:=-0.5+(LoopA-1)*1.0+0.5*CCPoints[LoopB];
      {Nudge Points}
      If (LoopB=1) Then
      Begin
        IntPoints[(LoopA-1)*PointNumber+LoopB]:=IntPoints[(LoopA-1)*PointNumber+LoopB]+1.0e-12;
      End Else If (LoopB=PointNumber) Then
      Begin
        IntPoints[(LoopA-1)*PointNumber+LoopB]:=IntPoints[(LoopA-1)*PointNumber+LoopB]-1.0e-12;
      End;
      {Map Weights}
      IntWeights[(LoopA-1)*PointNumber+LoopB]:=0.5*CCWeights[LoopB];
    End;
  End;

  {Set Results}
  SetLength(CCPoints,TotIntPoints+1);
  SetLength(CCWeights,TotIntPoints+1);
  For LoopA:=1 To TotIntPoints Do CCPoints[LoopA]:=IntPoints[LoopA];
  For LoopA:=1 To TotIntPoints Do CCWeights[LoopA]:=IntWeights[LoopA];
  PointNumber:=TotIntPoints;

End;


{Interpolation Of Table Data}
Function InterpolateTable(RowValue,ColValue: Double;RowCount,ColCount: Integer4;Table: Double2DArray): Double;
Var
  FoundRow,FoundCol: Boolean;
  IsLastRow,IsLastCol: Boolean;
  IsFirstRow,IsFirstCol: Boolean;
  FoundRowCount,FoundColCount: Integer4;
  XValues,YValues: Array2Double;
  TableValues: Array2x2Double;
  FirstColIndex,SecondColIndex: Integer4;
  FirstRowIndex,SecondRowIndex: Integer4;
  R1_Value,R2_Value: Double;
Begin
  {Find Row Values}
  FoundRow:=FALSE;
  FoundRowCount:=1;
  While (Not(FoundRow))And(FoundRowCount<RowCount) Do
  Begin
    Inc(FoundRowCount);
    FoundRow:=(Table[FoundRowCount,1]>RowValue);
  End;
  IsLastRow:=(FoundRowCount=RowCount)And(Not(FoundRow));
  IsFirstRow:=(FoundRowCount=2);
  {Find Columns Value}
  FoundCol:=FALSE;
  FoundColCount:=1;
  While (Not(FoundCol))And(FoundColCount<ColCount) Do
  Begin
    Inc(FoundColCount);
    FoundCol:=(Table[1,FoundColCount]>ColValue);
  End;
  IsLastCol:=(FoundColCount=ColCount)And(Not(FoundCol));
  IsFirstCol:=(FoundColCount=2);
  {Recover X Values}
  If (IsFirstCol)Or(IsLastCol) Then
  Begin
    FirstColIndex:=FoundColCount;
    SecondColIndex:=FoundColCount;
  End Else Begin
    FirstColIndex:=FoundColCount-1;
    SecondColIndex:=FoundColCount;
  End;
  {Recover Y Values}
  If (IsFirstRow)Or(IsLastRow) Then
  Begin
    FirstRowIndex:=FoundRowCount;
    SecondRowIndex:=FoundRowCount;
  End Else Begin
    FirstRowIndex:=FoundRowCount-1;
    SecondRowIndex:=FoundRowCount;
  End;
  {Store Values}
  XValues[1]:=Table[1,FirstColIndex];
  XValues[2]:=Table[1,SecondColIndex];
  YValues[1]:=Table[FirstRowIndex,1];
  YValues[2]:=Table[SecondRowIndex,1];
  {Recover Table Values}
  TableValues[1,1]:=Table[FirstRowIndex,FirstColIndex];
  TableValues[1,2]:=Table[SecondRowIndex,FirstColIndex];
  TableValues[2,1]:=Table[FirstRowIndex,SecondColIndex];
  TableValues[2,2]:=Table[SecondRowIndex,SecondColIndex];
  {Eval Interpolation}
  {Eval R1_Value}
  If (ABS(TableValues[1,1]-TableValues[2,1])<MathZero) Then R1_Value:=TableValues[1,1]
  Else R1_Value:=TableValues[1,1]*((XValues[2]-ColValue)/(XValues[2]-XValues[1]))+TableValues[2,1]*((ColValue-XValues[1])/(XValues[2]-XValues[1]));
  {Eval R2_Value}
  If (ABS(TableValues[1,2]-TableValues[2,2])<MathZero) Then R2_Value:=TableValues[1,2]
  Else R2_Value:=TableValues[1,2]*((XValues[2]-ColValue)/(XValues[2]-XValues[1]))+TableValues[2,2]*((ColValue-XValues[1])/(XValues[2]-XValues[1]));
  {Eval Final Result}
  If (ABS(R1_Value-R2_Value)<MathZero) Then Result:=R1_Value
  Else Result:=R1_Value*((YValues[2]-RowValue)/(YValues[2]-YValues[1]))+R2_Value*((RowValue-YValues[1])/(YValues[2]-YValues[1]));
End;

{Full Matrix Vector Multiplication}
Function FullMultiply_0Based(RowCount,ColCount: Integer4;
                             Mat: Double2DArray;
                             Vector: DoubleArray;
                             Var ResultantVector: DoubleArray): String;
Var
  LoopA,LoopB: Integer4;
Begin
  {Initialize Result}
  Result:='';
  {Check Dimension Of the Vector}
  If (Length(Vector)<>ColCount) Then
  Begin
    Result:='Error FullMultiply_0Based: Invalid Vector Dimensions.';
    Exit;
  End;
  For LoopA:=0 To (RowCount-1) Do
  Begin
    ResultantVector[LoopA]:=0.0;
    For LoopB:=0 To (ColCount-1) Do
    Begin
      ResultantVector[LoopA]:=ResultantVector[LoopA]+Mat[LoopA+1,LoopB+1]*Vector[LoopB];
    End;
  End;
End;

{Full Matrix Vector Multiplication}
Function FullMultiply_1Based(RowCount,ColCount: Integer4;
                             Mat: Double2DArray;
                             Vector: DoubleArray;
                             Var ResultantVector: DoubleArray): String;
Var
  LoopA,LoopB: Integer4;
Begin
  {Initialize Result}
  Result:='';

  {Allocate}
  SetLength(ResultantVector,RowCount+1);

  {Check Dimension Of the Vector}
  If (Length(Vector)<>ColCount+1) Then
  Begin
    Result:='Error FullMultiply_0Based: Invalid Vector Dimensions.';
    Exit;
  End;
  For LoopA:=1 To RowCount Do
  Begin
    ResultantVector[LoopA]:=0.0;
    For LoopB:=1 To ColCount Do
    Begin
      ResultantVector[LoopA]:=ResultantVector[LoopA]+Mat[LoopA,LoopB]*Vector[LoopB];
    End;
  End;
End;


{Full Matrix Vector Multiplication: Trasposed}
Function FullMultiply_T_0Based(RowCount,ColCount: Integer4;
                               Mat: Double2DArray;
                               Vector: DoubleArray;
                               Var ResultantVector: DoubleArray): String;
Var
  LoopA,LoopB: Integer4;
Begin
  {Initialize Result}
  Result:='';
  {Check Dimension Of the Vector}
  If (Length(Vector)<>RowCount) Then
  Begin
    Result:='Error FullMultiply_0Based: Invalid Vector Dimensions.';
    Exit;
  End;
  For LoopA:=0 To ColCount-1 Do
  Begin
    ResultantVector[LoopA]:=0.0;
    For LoopB:=0 To RowCount-1 Do
    Begin
      ResultantVector[LoopA]:=ResultantVector[LoopA]+Mat[LoopB+1,LoopA+1]*Vector[LoopB];
    End;
  End;
End;

{Full Matrix Vector Multiplication: Trasposed}
Function FullMultiply_T_1Based(RowCount,ColCount: Integer4;
                               Mat: Double2DArray;
                               Vector: DoubleArray;
                               Var ResultantVector: DoubleArray): String;
Var
  LoopA,LoopB: Integer4;
Begin
  {Initialize Result}
  Result:='';
  {Check Dimension Of the Vector}
  If (Length(Vector)<>RowCount+1) Then Result:='Error FullMultiply_0Based: Invalid Vector Dimensions.';
  For LoopA:=1 To ColCount Do
  Begin
    ResultantVector[LoopA]:=0.0;
    For LoopB:=1 To RowCount Do
    Begin
      ResultantVector[LoopA]:=ResultantVector[LoopA]+Mat[LoopB,LoopA]*Vector[LoopB];
    End;
  End;
End;

{Check Relative Difference Between two Vectors}
Function CheckSolutionDifferece(Size: Integer4;
                                FirstVec: DoubleArray;
                                SecondVec: DoubleArray): Double;
Var
  LoopA: Integer4;
  Norm: Double;
Begin
  Result:=0.0;
  Norm:=0.0;
  For LoopA:=1 To Size Do
  Begin
    Result:=Result+Sqr(FirstVec[LoopA]-SecondVec[LoopA]);
    Norm:=Norm+Sqr(FirstVec[LoopA]);
  End;
  Norm:=Sqrt(Norm);
  If (ABS(Norm)<MathZero) Then Result:=Sqrt(Result)
  Else Result:=Sqrt(Result)/Norm;
End;

{Normalize Columns}
Procedure TransformToUnitaryColums(RowCount,ColCount: Integer4;
                                   OrigAMat: Double2DArray;
                                   Var ColNorms: DoubleArray;
                                   Var Mat: Double2DArray);
Var
  LoopA,LoopB: Integer4;
Begin
  {Allocate}
  SetLength(Mat,RowCount+1);
  For LoopA:=1 To RowCount Do SetLength(Mat[LoopA],ColCount+1);
  SetLength(ColNorms,ColCount+1);

  {Eval Column Norm Vector}
  For LoopA:=1 To ColCount Do
  Begin
    ColNorms[LoopA]:=0.0;
    For LoopB:=1 To RowCount Do
    Begin
      ColNorms[LoopA]:=ColNorms[LoopA]+Sqr(OrigAMat[LoopB,LoopA]);
    End;
    ColNorms[LoopA]:=Sqrt(ColNorms[LoopA]);
  End;

  {Write Unitary Column Problem}
  For LoopA:=1 To RowCount Do
  Begin
    For LoopB:=1 To ColCount Do
    Begin
      If ABS(ColNorms[LoopB])>MathZero Then
      Begin
        Mat[LoopA,LoopB]:=OrigAMat[LoopA,LoopB]/ColNorms[LoopB];
      End Else Begin
        Mat[LoopA,LoopB]:=OrigAMat[LoopA,LoopB];
      End;
    End;
  End;
End;

{Eval Integration Points}
(*Function EvalGaussLegendreIntegratinoPoints(PointNumber: Integer4;Var CCPoints, CCWeights: DoubleArray): Boolean;
Begin
function I = gauss(f,n)             % (n+1)-pt Gauss quadrature of f

  For LoopA:=0 To PointNumber Do Beta[LoopA]=0.5/Sqrt(1.0-(1.0/Sqrt(2.0*LoopA)));

beta =  % 3-term recurrence coeffs
T = diag(beta,1) + diag(beta,-1);   % Jacobi matrix
[V,D] = eig(T);                     % eigenvalue decomposition
x = diag(D); [x,i] = sort(x);       % nodes (= Legendre points)
w = 2*V(1,i).^2;                    % weights
I = w*feval(f,x);                   % the integral
End;*)

Procedure PrintMatrixToCSV(MatFileName: String;TotalRows,TotalCols: Integer4;AMat: Double2DArray);
Var
  LoopA,LoopB: Integer4;
  MatFile: TextFile;
Begin
  AssignFile(MatFile,MatFileName);
  Rewrite(MatFile);
  For LoopA:=1 To TotalRows Do
  Begin
    For LoopB:=1 To TotalCols Do
    Begin
      Write(MatFile,'"'+FloatToStr(AMat[LoopA,LoopB])+'";');
    End;
    Writeln(MatFile);
  End;
  CloseFile(MatFile);
End;

{Check If An Integer Vector Is Null}
Function IsIntegerVectorNull(Size: Integer4;IntVector: Integer4Array): Boolean;
Var
  LoopA: Integer4;
Begin
  Result:=TRUE;
  For LoopA:=1 To Size Do
  Begin
    If (IntVector[LoopA]<>0) Then
    Begin
      Result:=FALSE;
      Exit;
    End;
  End;
End;

end.
