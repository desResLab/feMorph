#ifndef FEMFACE_H
#define FEMFACE_H

#include <stdio.h>
#include <vector>

#include "femNode.h"

class femFace{
public:
  // Data Members
  int number;
  int group;
  std::vector<int> faceNodes;
  std::vector<int> faceElements;
  // Constructor and Destructor
  femFace();
  femFace(std::vector<int> nodes);
  femFace(int tempNumber,std::vector<int> nodes);
  femFace(femFace* other);
  ~femFace();
  // MEMBER FUNCTIONS
  void evalFaceCentroid(std::vector<femNode*> nodeList, double* centroid);
  void evalFaceNormal(std::vector<femNode*> nodeList, double* centroid);
};

#endif // FEMFACE_H
