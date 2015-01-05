#ifndef FEMVECTOR_H
#define FEMVECTOR_H

# include "femModel.h"

class femVector{
  public:
    femDoubleVec values;

    // CONSTRUCTOR
    femVector(int total);

    // ASSEMBLE IN DENSE COLUMN FORMAT
    void assemble(femDoubleVec vec,femIntVec indices);
    void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
};

#endif // FEMVECTOR_H
