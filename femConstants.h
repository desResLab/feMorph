#ifndef FEMCONSTANTS_H
#define FEMCONSTANTS_H

// Constants: NUMERICAL
const double kPI = 3.141592653589793238462;

// Constants: NODES
const int kNodeDofs = 6;

// Constants: ELEMENTS
// Number of face in tetrahedra
const int kTetra4Nodes = 4;
// Number of face in tetrahedra
const int kTetraFaces = 4;

// Constants: ALGORITHMS
// Numerical Zero
const int kMathZero = 1.0e-7;
// Number Of Monte Carlo test for Starting Element
const int kMaxEnclosingMC = 100;

#endif // FEMCONSTANTS_H
