//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};

//+
Field[1] = Mean;
//+
MeshSize {2, 6, 8, 4, 1, 3, 5, 7} = 0.01;
//+
Delete Field [1];
//+
MeshSize {3, 2, 7, 1, 5, 8, 6, 4} = 0.05;
//+
MeshSize {7, 5, 3, 8, 1, 6, 4, 2} = 0.05;
