#ifndef FEMTYPES_H
#define FEMTYPES_H

#include <vector>

using namespace std;

// TYPE OF RESULT
enum femResultType{frNode,frElement};

// ENTITY DIMENSION
enum elDim{d1,d2,d3};

// TYPE OF INTEGRATION RULES
enum intRuleType{irFirstOrder,irSecondOrder};

typedef vector<vector<double> > femDoubleMat;
typedef vector<vector<vector<double> > > femDoubleDOFMat;
typedef vector<vector<double> > femDoubleMat;
typedef vector<vector<int> > femIntMat;
typedef vector<double> femDoubleVec;
typedef vector<int> femIntVec;

#endif // FEMTYPES_H
