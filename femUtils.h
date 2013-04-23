#ifndef FEMUTILS_H
#define FEMUTILS_H

// Random Number Generator
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace femUtils{

// GLOBAL VARIABLE FOR GENERATION
static boost::random::mt19937 intGen;

// Generate Uniform Integer (global Variable Defined)
inline int GenerateUniformIntegers(int lowIdx, int upIdx) {
    boost::random::uniform_int_distribution<> dist(lowIdx, upIdx);
    return dist(intGen);
}

// Extral product between vectors
inline void Do3DExternalProduct(double* v1,double* v2,double* v3){
  v3[0] = v1[1]*v2[2]-v2[1]*v1[2];
  v3[1] = v1[2]*v2[0]-v2[2]*v1[0];
  v3[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

// Eval Internal Product
inline double Do3DInternalProduct(double* v1, double* v2){
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

// ===================
// Normalize 3D Vector
// ===================
inline void Normalize3DVector(double* vector){
  // Define modulus
  double modulus = sqrt((vector[0]*vector[0])+(vector[1]*vector[1])+(vector[2]*vector[2]));
  // Normalize
  if (modulus>0.0){
    vector[1] = vector[1]/modulus;
  }
  if (modulus>0.0){
    vector[2] = vector[2]/modulus;
  }
  if (modulus>0.0) {
    vector[3] = vector[3]/modulus;
  }
}

// =================================================
// General Rotation Around a Vector with quaternions
// The Angle is in degrees
// =================================================
inline void Rotate3DVectorAroundAxis(double* vec, double angle,double* axis){
  double quatMat[4][4];
  double quatVect[4] = {0.0};
  double quatRes[4] = {0.0};
  // Normalize Axis
  Normalize3DVector(axis);
  // Write the compoents of the quaternion
  double w = cos(0.5*angle*(kPI/180.0));
  double x = axis[0]*sin(0.5*angle*(kPI/180.0));
  double y = axis[1]*sin(0.5*angle*(kPI/180.0));
  double z = axis[2]*sin(0.5*angle*(kPI/180.0));
  // Eval rotation Matrix
  quatMat[0][0] = w*w+x*x-y*y-z*z;
  quatMat[0][1] = 2.0*x*y+2.0*w*z;
  quatMat[0][2] = 2.0*x*z-2.0*w*y;
  quatMat[0][3] = 0.0;
  quatMat[1][0] = 2.0*x*y-2.0*w*z;
  quatMat[1][1] = w*w-x*x+y*y-z*z;
  quatMat[1][2] = 2.0*y*z+2.0*w*x;
  quatMat[1][3] = 0.0;
  quatMat[2][0] = 2.0*x*z+2.0*w*y;
  quatMat[2][1] = 2.0*y*z-2.0*w*x;
  quatMat[2][2] = w*w-x*x-y*y+z*z;
  quatMat[2][3] = 0.0;
  quatMat[3][0] = 0.0;
  quatMat[3][1] = 0.0;
  quatMat[3][2] = 0.0;
  quatMat[3][3] = w*w+x*x+y*y+z*z;
  // Expand vector to rotate
  quatVect[3] = 0.0;
  for(int loopA=0;loopA<3;loopA++){
    quatVect[loopA] = vec[loopA];
  }
  // Matrix vector multiplication
  for(int loopA=0;loopA<4;loopA++){
    quatRes[loopA] = 0.0;
    for(int loopB=0;loopB<4;loopB++){
      quatRes[loopA] += quatVect[loopB] * quatMat[loopB][loopA];
    }
  }
  // Final Copy of the Result
  for(int loopA=0;loopA<3;loopA++){
    vec[loopA] = quatRes[loopA];
  }
}


}

#endif // FEMUTILS_H
