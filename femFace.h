#ifndef FEMFACE_H
#define FEMFACE_H

#include <stdio.h>
#include <vector>

class femFace
{
public:
  // Data Members
  std::vector<int> faceNodes;
  std::vector<int> faceElements;
  // Constructor and Destructor
  femFace();
  femFace(std::vector<int> nodes);
  femFace(femFace* other);
  ~femFace();
};

#endif // FEMFACE_H