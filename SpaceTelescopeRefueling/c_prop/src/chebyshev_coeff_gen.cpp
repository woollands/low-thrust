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
#include <complex.h>
#include "const.hpp"
#include "c_functions.hpp"
#include "lsq_chebyshev_fit.h"
#include "../../../../../../Dropbox/CODE/cspice/src/cspice/spkezr_c.c"

void chebyshev_coeff_gen( int N, double t0, double tf, double coeff ){

  // Chebyshev Least Squares Matrices
  int s = -1.0;
  // Least Squares Operator (A)
  double Ta[N*N];
  memset( Ta, 0.0, (N*N*sizeof(double)));
  double A[N*N];
  memset( A, 0.0, (N*N*sizeof(double)));
  lsq_chebyshev_fit(s,N,N,Ta,A);


  // Compute Chebyshev Cosine Nodes
  double w1 = (tf+t0)/2.0;
  double w2 = (tf-t0)/2.0;
  double tau[N+1];
  memset( tau, 0.0, ((N+1)*sizeof(double)));
  double time[(N+1)];
  memset( time, 0.0, ((N+1)*sizeof(double)));
  for (int i=0; i<=N; i++){
      tau[i]   = -cos(i*C_PI/N);
      times[i] = tau[i]*w2 + w1;
    }

  // Compute States
  double states[6] = {0.0};
  double owlt = 0.0;
  for (i=0; i<=N; i++){
    spkezr_c("SUN", times[i], "ECLIPJ2000", "LT+S", "EARTH", states, &owlt);
    for (int j=0; j<=2; j++){
      X[ID2(i-1,j-1,N+1)] = states[j];
    }
  }

  // Compute Chebyshev Coefficients
  matmul(coeff,A,X,N+1,N+1,3,N+1,N+1,N+1);

  return;

}
