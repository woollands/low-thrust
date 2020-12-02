/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     MDec 2020
*  LAST MODIFIED:    Dec 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Compute Chebyshev coefficients from cspice states
*  REFERENCE:
*/

// #include "../../../../cspice/src/cspice/conics_c.c"
# include "../../../../../../Dropbox/CODE/cspice/src/cspice/furnsh_c.c"
// # include "../../../../../../Dropbox/CODE/cspice/src/cspice/str2et_c.c"

#include "const.hpp"
#include "c_functions.hpp"
#include "chebyshev_coeff_gen.hpp"
#include "ephem_gen_chebyshev.hpp"

int main()
{

  // Load Spice Kernels
  furnsh_c ( "../../../../../Dropbox/CODE/mice/kernels/naif0012.tls" );
  furnsh_c ( "../../../../../Dropbox/CODE/mice/kernels/de438.bsp" );

  double t0 = 0.0;
  double tf = 100.0*86400.0;

  //////////////////////////////////////////////////////////////////////////////
  // Generate Chebyshev coefficients over specified time range (run once at the beginning of code)
  int N;
  N = (int) ((tf-t0)/30/86400) * 40;
  if (N == 0){
    N = 40;
  }
  double coeff[(N+1)*6];
  memset( coeff, 0.0, ((N+1)*6*sizeof(double)));
  chebyshev_coeff_gen(N,t0,tf,coeff);
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // Compute states using previously computed Chebyshev coefficients
  double t = 17.0*86400.0;
  double states[6] = {0.0};
  ephem_gen_chebyshev(N,t0,tf,coeff,t,states); // Call this from eclipse function

  printf("Robyn's version \n");
  for (int i=0; i<=5; i++){
    printf("%f\t",states[i]);
  }
  printf("\n");

  //////////////////////////////////////////////////////////////////////////////
  // Check using cspice
  double states1[6] = {0.0};
  double owlt = 0.0;
  spkezr_c("SUN", t, "ECLIPJ2000", "LT+S", "EARTH", states1, &owlt);

  printf("CSPICE version \n");
  for (int i=0; i<=5; i++){
    printf("%f\t",states1[i]);
  }
  printf("\n");

}
