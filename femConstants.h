#ifndef FEMCONSTANTS_H
#define FEMCONSTANTS_H

// Constants: NUMERICAL
const double kPI = 3.141592653589793238462;

// Constants: NODES
const int kNodeDofs = 6;

// Constants: Internal Grid
const int kGridSizeX = 1;
const int kGridSizeY = 1;
const int kGridSizeZ = 1;

// Constants: ELEMENTS
// Number of face in tetrahedra 4
const int kTetra4Nodes = 4;
// Number of face in tetrahedra 10
const int kTetra10Nodes = 10;
// Number of face in tetrahedra
const int kTetraFaces = 4;

// Constants: ALGORITHMS
// Numerical Zero
const double kMathZero = 1.0e-7;
// Number Of Monte Carlo test for Starting Element
const int kMaxEnclosingMC = 5;
// Grid Tolerance
const double kGridTol = 1.0e-3;


#endif // FEMCONSTANTS_H
