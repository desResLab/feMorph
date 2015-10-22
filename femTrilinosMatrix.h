#ifndef FEMTRILINOSMATRIX_H
#define FEMTRILINOSMATRIX_H

# include "femMatrix.h"

# include "Epetra_Map.h"
# include "Epetra_BlockMap.h"
# include "Epetra_FEVbrMatrix.h"

# include "femModel.h"
# include "femException.h"

class femTrilinosMatrix : public femMatrix{
  private:
    Epetra_Map* Map;
    Epetra_BlockMap* BlockMap;
    // Graph used to initialize the topology
    Epetra_CrsGraph* MatGraph;
    // Block Finite element matrix in Trilinos
    Epetra_FEVbrMatrix* MatValues;
  public:
    femTrilinosMatrix(femModel* model);

    // VIRTUAL FUNCTIONS
    virtual void assemble(femDoubleMat elMat,femIntVec connections);
    virtual void applyDirichelet(femIntVec dofs);
    // I/O
    virtual void writeToFile(string fileName);
    virtual double getRowSum(int loopA);
    virtual void clearRowAndColumn(int dof);

    // EXTRA
    void assembleDOF(femDoubleDOFMat elMat,femIntVec connections);

};

#endif // FEMTRILINOSMATRIX_H
