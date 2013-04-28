#ifndef FEMPOINT_H
#define FEMPOINT_H

class femPoint
{
public:
  // Data Members
  double coords[3];
  // Constructor and Distructor
  femPoint(double* localCoords);
  ~femPoint();
};

#endif // FEMPOINT_H
