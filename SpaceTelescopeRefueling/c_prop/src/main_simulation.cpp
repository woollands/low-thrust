/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Mar 2020
*  LAST MODIFIED:    Mar 2020
*  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
*  DESCRIPTION:      Integrate low thrust trajectory using this the BOOST library's bulirsch-stoer integrator
*  REFERENCE:        Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE_1_0.txt
*                    or copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

#include <boost/array.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/config.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

#include <boost/math/tools/roots.hpp>

// #include "../../../../cspice/src/cspice/conics_c.c"
# include "../../../../../../Dropbox/CODE/cspice/src/cspice/furnsh_c.c"
# include "../../../../../../Dropbox/CODE/cspice/src/cspice/str2et_c.c"

#include "const.hpp"
#include "states_twobody_thrust.hpp"
#include "chebyshev_coeff_gen.hpp"

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array< double , 14 > state_type;

ofstream out;

void write_out( const state_type &x , const double t )
{
    out << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << '\t' << x[4] << '\t' << x[5] << '\t' << x[6] <<  '\t' << x[7] << '\t' << x[8] << '\t' << x[9] << '\t' << x[10] << '\t' << x[11] << '\t' << x[12] << '\t' << x[13] << '\t' << endl;

    #include "spacecraft_params.hpp"
    #include "eclipse_model.hpp"

    double Pa;
    double TF;
    double y[7] = {0.0};
    for (int i=0; i<=6; i++){
      y[i] = x[i];
    }
    eclipse_model(t,y,P,true,&Pa,&TF);
    // printf("TF %f\n",TF);

}

double rho = 1.0;

int main()
{

    // Load Spice Kernels
    furnsh_c ( "../../../../../Dropbox/CODE/mice/kernels/naif0012.tls" );
    furnsh_c ( "../../../../../Dropbox/CODE/mice/kernels/de438.bsp" );

    bulirsch_stoer< state_type > controlled_stepper( 1.0E-8 , 1.0E-14 );

    state_type x = {{0.0}};

    for (int i=0; i<1; i++){

      // rho = 0.1*rho;
      // Initial Conditions
      x[0] = 1.82260259877705;
      x[1] = 0.725;
      x[2] = 0.0;
      x[3] = 0.0611626201504845;
      x[4] = 0.0;
      x[5] = 0.0;
      x[6] = 100.0;
      x[7] = -4.75169058027082;
      x[8] = -12.6007096436135;
      x[9] = 0.214505265997245;
      x[10] = 5.95508267111558;
      x[11] = -0.0366858149807001;
      x[12] = 0.00305135494989929;
      x[13] = 0.118203061064424;

      // Time (Canonical Units)
      double t0 = 0.0;
      double tf = 642.549910873478;
      double dt = 1.0;

      // out.open( "./output/data_extra.dat" );
      out.open( "./output/data.dat" );
      out.precision(16);
      // integrate_adaptive( controlled_stepper , states_twobody_thrust , x , t0 , tf , dt , write_out); // Use for accurate final time
      integrate_const( controlled_stepper , states_twobody_thrust, x , t0 , tf , dt , write_out ); // Use for plotting
      out.close();

    }

    // ConstSpiceChar * str;
    // SpiceDouble    * et;
    // str = "1 JAN 2020 UTC";
    // str2et_c(str,et);

    double t0 = 0.0;
    double tf = 642.549910873478;
    int N;
    N = (int) ((tf-t0)/30/86400) * 40;
    if (N == 0){
      N = 40;
    }
    double coeff[N*N];
    memset( coeff, 0.0, (N*N*sizeof(double)));
    chebyshev_coeff_gen(N,t0,tf,coeff);
    // spkezr_c("SUN", t, "J2000", "LT+S", "EARTH", states, &owlt);

}
