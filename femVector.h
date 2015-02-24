#ifndef FEMVECTOR_H
#define FEMVECTOR_H

# include "femModel.h"

class femVector{
  public:
    femDoubleVec values;

    // CONSTRUCTOR
    femVector(int total);

    // GETTER AND SETTER
    int getSize(){return (int)values.size();}

    // ASSEMBLE IN DENSE COLUMN FORMAT
    void assemble(femDoubleVec vec,femIntVec indices);
    void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
    void writeToFile(string fileName);
};

#endif // FEMVECTOR_H
