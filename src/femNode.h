#ifndef FEMNODE_H
#define FEMNODE_H

class femNode
{
public:
    // Data Members
    int nodeNumber;
    double coords[3];
    double displacements[6];
    double force[3];
    // Constructor and Destructor
    femNode(int number, double coordX, double coordY, double coordZ);
    femNode(int number, double* coords, double* disps);
    ~femNode();
    // Getter and Setter Functions
    double* getNodeDisplacements(){return displacements;}
    // Member Functions
    void setDisplacements(double DX, double DY, double DZ, double RX, double RY, double RZ);
    // Node Transformation
    void TransformNodeCoords(double* origin, double** rotMat, double* newCoords, const int transformType, const int dispType, double dispFactor);
    void TransformDisplacements(double** rotMat, double* newDisps);
};

#endif // FEMNODE_H
