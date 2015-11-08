#ifndef FEMTRILINOSMATRIX_H
#define FEMTRILINOSMATRIX_H

# include <iostream>
# include <fstream>

# include "../femMatrix.h"

#ifdef USE_MPI
  # include "mpi.h"
  # include "Epetra_MpiComm.h"
#else
  # include "Epetra_SerialComm.h"
#endif
# include "Epetra_Map.h"
# include "Epetra_BlockMap.h"
# include "Epetra_FEVbrMatrix.h"

# include "../femModel.h"
# include "../femException.h"

class femTrilinosMatrix : public femMatrix{
  public:
    #ifdef USE_MPI
      Epetra_MpiComm* EpetraComm;
    #else
      Epetra_SerialComm* EpetraComm;
    #endif

    // Epetra_FEVector
    int nodeDOFs;
    Epetra_BlockMap* map;
    Epetra_CrsGraph* graph;
    Epetra_FEVbrMatrix* values;

    // Construct the matrix
    femTrilinosMatrix(femModel* model,int nodeDOFs);
    virtual ~femTrilinosMatrix();

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat nodeMat,femIntVec elConnections);
    virtual void blockAssemble(femDoubleBlockMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    void         applyBlockDirichelet(femIntVec gNodesIdx, int dof);
    // I/O
    virtual void   writeToFile(string fileName);
    virtual double getRowSum(int loopA);
    virtual void   clearRowAndColumn(int dof);

    // COMPLETE FILLING: TRILINOS ONLY MATRIX
    virtual void completeFill();
};

#endif // FEMTRILINOSMATRIX_H
