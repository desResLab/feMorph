#ifndef FEMCONSTANTS_H
#define FEMCONSTANTS_H

// Constants: NUMERICAL
const double kPI = 3.141592653589793238462;

// Constants: NUMBER OF DIMENSIONS
const int kDims = 3;

// Constants: NODES
const int kNodeDofs = 6;

// Constants: transformation
const int kDirect = 0;
const int kReverse = 1;
// Constants: Add displacements
const int kUndeformed = 0;
const int kDeformed = 1;

// Constants: Internal Grid
const int kGridSizeX = 10;
const int kGridSizeY = 10;
const int kGridSizeZ = 10;

// Constants: ELEMENTS
// Number of face in tetrahedra 4
const int kTri3Nodes = 3;
// Number of face in tetrahedra 4
const int kTetra4Nodes = 4;
// Number of face in tetrahedra 10
const int kTetra10Nodes = 10;
// Number of face in tetrahedra
const int kTetraFaces = 4;

// Constants: Model Quality
// Volume
const int kVolume = 0;
// Mixed Product
const int kMixedProduct = 1;


// Constants: ALGORITHMS
// Numerical Zero
const double kMathZero = 1.0e-7;
// Number Of Monte Carlo test for Starting Element
const int kMaxEnclosingMC = 5;
// Grid Tolerance
const double kGridTol = 1.0e-3;
// Number of slices in stenosis
const int kStenosisSlices = 20;
// Tolerance in stenotic value
const double kStenosisTolerance = 1.0e-4;
// Max Picard Iterations for Stenosis
const int kMaxPicardStenosisIt = 20;
// Transverse Scaling factor
const double kTransverseScalingFactor = 1.05;
// Determine nodes to inclue in stenosis box
const double kStenosisBoxFactor = 0.45;
// Angle Limit for normal compatibility check
const double kNormalAngleLimit = 50.0;

#endif // FEMCONSTANTS_H
