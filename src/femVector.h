#ifndef FEMVECTOR_H
#define FEMVECTOR_H

# include "femModel.h"

# include "femException.h"

class femVector{
  public:
    // CONSTRUCTOR
    femVector();

    // GETTER AND SETTER
    // Get total size of the array
    virtual int getSize();
    // Get single component
    virtual double getComponent(int id);
    virtual void setComponent(int id, double entry);
    virtual void addComponent(int id, double entry);

    // ASSEMBLE IN DENSE COLUMN FORMAT
    virtual void assemble(femDoubleVec vec,femIntVec indices);
    virtual void blockAssemble(femDoubleBlockVec vec,femIntVec indices);
    virtual void GlobalAssemble();
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
    virtual double getComponent(int id);
    virtual void setComponent(int id, double entry);
    virtual void addComponent(int id, double entry);

    // ASSEMBLE IN DENSE COLUMN FORMAT
    virtual void assemble(femDoubleVec vec,femIntVec indices);
    virtual void blockAssemble(femDoubleBlockVec vec,femIntVec indices);
    virtual void GlobalAssemble();
    virtual void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
    virtual void writeToFile(string fileName);
};

#endif // FEMVECTOR_H
