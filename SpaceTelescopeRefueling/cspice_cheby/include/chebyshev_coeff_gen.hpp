/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Dec 2020
*  LAST MODIFIED:    Dec 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Header file
*/

#ifndef __CCG__
#define __CCG__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.hpp"

void chebyshev_coeff_gen( int N, double t0, double tf, double* coeff);

#endif
