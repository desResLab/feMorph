# include "femTrilinosMatrix.h"

// CONSTRUCTOR
femTrilinosMatrix::femTrilinosMatrix(femModel* model){
    throw femException("Not Implemented.\n");
  // Create Matrix Topology From Model

    // Epetra_CrsGraph::Epetra_CrsGraph	(	Epetra_DataAccess 	CV,
    // const Epetra_BlockMap & 	RowMap,
    // const int * 	NumIndicesPerRow,
    // bool 	StaticProfile = false

}


// VIRTUAL FUNCTIONS
void  femTrilinosMatrix::assemble(femDoubleMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}
void  femTrilinosMatrix::assembleDOF(femDoubleDOFMat elMat,femIntVec connections){
  throw femException("Not Implemented.\n");
}
void  femTrilinosMatrix::applyDirichelet(femIntVec dofs){
  throw femException("Not Implemented.\n");
}
// I/O
void   femTrilinosMatrix::writeToFile(string fileName){
  throw femException("Not Implemented.\n");
}
double femTrilinosMatrix::getRowSum(int loopA){
  throw femException("Not Implemented.\n");
}
void   femTrilinosMatrix::clearRowAndColumn(int dof){
  throw femException("Not Implemented.\n");
}
