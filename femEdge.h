#ifndef FEMEDGE_H
#define FEMEDGE_H

#include <vector>

class femEdge
{
public:
  // Data Members
  int number;
  std::vector<int> edgeNodes;
  std::vector<int> edgeElements;
  // Constructor and Distructor
  femEdge();
  femEdge(int tempNumber,std::vector<int> nodes);
};

#endif // FEMEDGE_H
