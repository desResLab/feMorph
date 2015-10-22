#ifndef FEMVECTOR_H
#define FEMVECTOR_H

# include "femModel.h"

# include "femException.h"


class femVector{
  public:
    // CONSTRUCTOR
    femVector();

    // GETTER AND SETTER
    virtual int getSize(){return 0;}

    // ASSEMBLE IN DENSE COLUMN FORMAT
    virtual void assemble(femDoubleVec vec,femIntVec indices);
    virtual void assembleDOF(femDoubleDOFVec vec,femIntVec indices);
    virtual void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
    virtual void writeToFile(string fileName);
};

class femDenseVector: public femVector{
  public:
    femDoubleVec values;

    // CONSTRUCTOR
    femDenseVector(int total);

    // GETTER AND SETTER
    virtual int getSize(){return (int)values.size();}

    // ASSEMBLE IN DENSE COLUMN FORMAT
    virtual void assemble(femDoubleVec vec,femIntVec indices);
    virtual void assembleDOF(femDoubleDOFVec vec,femIntVec indices);
    virtual void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
    virtual void writeToFile(string fileName);
};

#endif // FEMVECTOR_H
