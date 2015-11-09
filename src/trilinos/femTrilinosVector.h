#ifndef FEMTRILINOSVECTOR_H
#define FEMTRILINOSVECTOR_H

# include <iostream>
# include <fstream>

# include "Epetra_FEVector.h"
# include "../femVector.h"
#ifdef USE_MPI
  # include "mpi.h"
  # include "Epetra_MpiComm.h"
#else
  # include "Epetra_SerialComm.h"
#endif

class femTrilinosVector : public femVector{
public:
  #ifdef USE_MPI
    Epetra_MpiComm* EpetraComm;
  #else
    Epetra_SerialComm* EpetraComm;
  #endif

  // Epetra_FEVector
  int nodeDOFs;
  Epetra_BlockMap* map;
  Epetra_FEVector* values;

  // CONSTRUCTOR
  femTrilinosVector(int totNodesInProc,int* localToGlobal,int nodeDOFs);
  femTrilinosVector(Epetra_BlockMap map,int totComponents,int nodeDOFs);

  // GETTER AND SETTER
  virtual int getSize();
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

#endif // FEMTRILINOSVECTOR_H
