#ifndef FEMTRILINOSVECTOR_H
#define FEMTRILINOSVECTOR_H

# include "femVector.h"

class femTrilinosVector : public femVector{
public:
  femDoubleDOFVec values;

  // CONSTRUCTOR
  femTrilinosVector(int total,int nodeDOFs);

  // GETTER AND SETTER
  virtual int getSize();

  // ASSEMBLE IN DENSE COLUMN FORMAT
  virtual void assemble(femDoubleVec vec,femIntVec indices);
  virtual void assembleDOF(femDoubleDOFVec vec,femIntVec indices);
  virtual void applyDirichelet(femIntVec diricheletBCNode,femDoubleVec diricheletBCValues);
  virtual void writeToFile(string fileName);

};

#endif // FEMTRILINOSVECTOR_H
