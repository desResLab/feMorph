#ifndef FEMTYPES_H
#define FEMTYPES_H

#include <vector>

// TYPE OF RESULT
enum femResultType{frNode,frElement};

// TYPE OF INTEGRATION RULES
enum intRuleType{irFirstOrder,irSecondOrder};

typedef std::vector<std::vector<double>> femDoubleMat;
typedef std::vector<std::vector<int>> femIntMat;
typedef std::vector<double> femDoubleVec;
typedef std::vector<int> femIntVec;

#endif // FEMTYPES_H
