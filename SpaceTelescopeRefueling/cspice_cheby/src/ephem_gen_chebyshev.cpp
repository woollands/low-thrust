/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com), David Gardiner (davegard@vt.edu)
*  DATE WRITTEN:     Dec 2020
*  LAST MODIFIED:    Dec 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Generates Chebyshev coefficients from SPICE data
*  REFERENCE:
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include "const.hpp"
#include "c_functions.hpp"

void ephem_gen_chebyshev( int N, double t0, double tf, double* coeff, double t, double* states){

  // Compute tau vector
  double w1 = (tf+t0)/2.0;
  double w2 = (tf-t0)/2.0;
  double tau = (t - w1)/w2;

  // Interpolate Chebyshev Polynomials
  double T[(N+1)];
  memset( T, 0.0, ((N+1)*sizeof(double)));
  for (int i=0; i<=N; i++){
    T[ID2(1,i+1,1)] = cos(i*acos(tau));
    // printf("T %f\n",T[ID2(1,i+1,1)]);
  }

  // Compute states
  matmul(T,coeff,states,1,N+1,6,1,N+1,1);

  return;

}
