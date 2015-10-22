#include "femEdge.h"

// Standard Constructor
femEdge::femEdge(){

}

// Overloaded Constructor
femEdge::femEdge(int tempNumber,std::vector<int> nodes){
  // Assign number
  number = tempNumber;
  // Assign Nodes
  for(unsigned int loopA=0;loopA<nodes.size();loopA++){
    edgeNodes.push_back(nodes[loopA]);
  }
}
