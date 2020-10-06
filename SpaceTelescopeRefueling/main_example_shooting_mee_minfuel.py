## @package main
#  \brief     Simple shooting method.
#  \details   Solves a simple problem to demonstrate the shooting method.
#  \author    Robyn Woollands
#  \version   1.0
#  \date      2020-09-17
#  \pre       Requires SPICE Toolkit (https://naif.jpl.nasa.gov/naif/toolkit.html)
#  \bug       No bugs known
#  \warning   May fail to converge for poor initial costate guess
#  \copyright Robyn Woollands

import pandas as pd
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import spiceypy as sp
from scipy.integrate import solve_ivp as integrator
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from matplotlib import rc
import const as cn
import sys
from mps import mps_mee_ocp
from eom import eom_mee_twobody_minfuel
from functions import classical2mee
from functions import mee2rv
from functions import switch_function
from plot_mee_minfuel import plot_mee_minfuel

# Load Spice Kernels
sp.furnsh( "/Users/robynmw/Dropbox/CODE/mice/kernels/naif0012.tls" )
sp.furnsh( "/Users/robynmw/Dropbox/CODE/mice/kernels/de438.bsp" )

def residuals(p0,tspan,states0,statesf,rho):
    ## Computes the error between the current and desired final state.
    states_in = np.zeros(14)
    states_in[0:7] = states0[0:7]
    states_in[7:14] = p0[0:7]
    # sol = integrator(eom_mee_twobody_minfuel,(tspan[0],tspan[-1]),states_in,method='LSODA',rtol=1e-12)
    sol = integrator(lambda t,y: eom_mee_twobody_minfuel(t,y,rho),tspan,states_in,method='LSODA',rtol=1e-12)
    # Free final mass and free final true longitude (angle)
    res = np.linalg.norm(statesf[0:5] - sol.y[0:5,-1])
    return res

# Time
tspan = np.array([0,642.54]) # Nondimenional Units
# Initial State
sma = 24505/cn.DU
ecc = 0.725
inc = 7
Om  = 0
w   = 0
M   = 0
m0  = 100
[p,f,g,h,k,L] = classical2mee(sma,ecc,np.deg2rad(inc),Om,w,M)
mee0 = [p, f, g, h, k, L, m0]
# Final State
sma = 42165/cn.DU
ecc = 0
inc = 0
Om  = 0
w   = 0
M   = 0
[p,f,g,h,k,L] = classical2mee(sma,ecc,np.deg2rad(inc),Om,w,M)
meef = [p, f, g, h, k, L]
# Costate Guess
p0 = np.zeros(7)
p0[0] = -4.620949386264961
p0[1] = -12.037907266888872
p0[2] = 0.208961742408408
p0[3] = 5.946490544006020
p0[4] = 0.042440374825155
p0[5] = 0.002470486378745
p0[6] = 0.117893582793439
# p0 = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01])

# Solve Min-Fuel Low Thrust
rho = 1 # Continuation Parameter (engine throttle)
p0_guess = p0
optres = minimize(residuals,p0_guess[0:7],args=(tspan,mee0,meef,rho),method='Nelder-Mead',tol=1e-2)
print("opt_err:",optres.fun)
# print(states0_new)
p0 = optres.x

sys.exit()
# Continuation on rho
rho      = 1
rho_rate = 0.9
mps_tol  = 1e-4
iter     = 0
while rho > 1e-4:

    if rho < 1e-4:
        mps_tol  = 1e-5
        rho_rate = 0.95
    # Method of Particular Solutions
    [soln,p0,mps_err] = mps_mee_ocp(tspan,p0,mee0,meef,rho,mps_tol)
    print("rho:",rho,"mps err:",mps_err)
    rho  = rho*rho_rate

    # Update
    iter = iter + 1

# Final Solution
rho = rho/rho_rate # rho for converged solution
mee0_new = np.zeros(14)
mee0_new[0:7] = mee0[0:7]
mee0_new[7:14] = p0[0:7]

# Check Solution
sol   = integrator(lambda t,y: eom_mee_twobody_minfuel(t,y,rho),tspan,mee0_new,method='LSODA',rtol=1e-13)
data  = np.array(soln.y,dtype=float).T
time  = np.array(soln.t,dtype=float).T
plot_mee_minfuel(time,data,rho)
