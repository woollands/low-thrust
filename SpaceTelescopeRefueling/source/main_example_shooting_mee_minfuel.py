## @package main_example_shooting_mee_minfuel
#  \brief     Main file for mee, low-thrust, min-fuel, orbit transfer problem setup
#  \details   Example shooting method using scipy.optimize and using the
#  method of particular solutions.
#  \author    Robyn Woollands
#  \pre       numpy, const.py
#  \bug       No bugs known

import pandas as pd
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import spiceypy as sp
from scipy.integrate import solve_ivp as integrator
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import const as cn
import sys
from mps import mps_mee_ocp
from eom import eom_mee_twobodyJ2_minfuel
from functions import classical2mee
from functions import mee2rv
from functions import switch_function
from plot_routines import plot_mee_minfuel
from mee_ocp_finite_difference import mee_ocp_central_difference
import spacecraft_params as sc

# Load Spice Kernels
sp.furnsh( "/Users/robynmw/Dropbox/CODE/mice/kernels/naif0012.tls" )
sp.furnsh( "/Users/robynmw/Dropbox/CODE/mice/kernels/de438.bsp" )

def mee_ocp_central_difference(x):
    """
    Central different computation of the Jacobian of the cost funtion.
    Parameters:
    ===========
    tspan   -- vector containting initial and final time
    p0      -- initial costates guess
    states0 -- initial state vector (modified equinoctial elements)
    rho     -- switch smoothing parameter
    eclipse -- boolean (true or false)
    Returns:
    ========
    Jac -- Jacobian
    External:
    =========
    numpy, spipy.integrate
    """

    # TEMP
    tspan = np.array([0,642.54])
    rho = 1
    eclipse = False
    meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])

    # Initialization
    IC1  = np.zeros(14)
    IC2  = np.zeros(14)
    grad = np.zeros(7)
    DEL  = 1e-5*(x[7:14]/np.linalg.norm(x[7:14]))

    for i in range(0,7):

        IC1 = x
        IC2 = x
        IC1[i+7] = IC1[i+7]+DEL[i]
        IC2[i+7] = IC2[i+7]-DEL[i]

        # Integrate Dynamics
        sol1 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC1,method='LSODA',rtol=1e-13)
        sol2 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC2,method='LSODA',rtol=1e-13)

        res1 = np.linalg.norm(meef[0:5] - sol1.y[0:5,-1])
        res2 = np.linalg.norm(meef[0:5] - sol2.y[0:5,-1])

        grad[i] = (res1-res2)/(2*DEL[i])

    return grad

def residuals_mee_ocp(p0,tspan,mee0,meef,rho,eclipse):
    """
    Computes error between current final state and desired final state for the
    low thrust, min-fuel, optimal trajectory using modified equinoctial elements.
    Parameters:
    ===========
    p0    -- current initial costate
    tspan -- vector containing initial and final time
    mee0  -- initial states (modified equinoctial elements)
    meef  -- final states (modified equinoctial elements)
    rho   -- switch smoothing parameter
    Returns:
    ========
    res -- residual
    External:
    =========
    numpy, scipy.integrate
    """
    ## Computes the error between the current and desired final state.
    states_in = np.zeros(14)
    states_in[0:7] = mee0[0:7]
    states_in[7:14] = p0[0:7]
    # sol = integrator(eom_mee_twobody_minfuel,(tspan[0],tspan[-1]),states_in,method='LSODA',rtol=1e-12)
    sol = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,states_in,method='LSODA',rtol=1e-12)
    # Free final mass and free final true longitude (angle)
    res = np.linalg.norm(meef[0:5] - sol.y[0:5,-1])
    print(res)
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
# p0[0] = -4.620949386264961
# p0[1] = -12.037907266888872
# p0[2] = 0.208961742408408
# p0[3] = 5.946490544006020
# p0[4] = 0.042440374825155
# p0[5] = 0.002470486378745
# p0[6] = 0.117893582793439
p0 = 0.1*np.random.rand(7)

# Solve Min-Fuel Low Thrust
rho = 1 # Continuation Parameter (engine throttle)
eclipse = False
p0_guess = p0
# optres = minimize(residuals_mee_ocp,p0_guess[0:7],args=(tspan,mee0,meef,rho),method='Nelder-Mead',tol=1e-2)
optres = minimize(residuals_mee_ocp,p0_guess[0:7],args=(tspan,mee0,meef,rho,eclipse),method='BFGS',tol=1e-2)
print("opt_err:",optres.fun)
# print(states0_new)
p0 = optres.x

sys.exit()
# Continuation on rho
rho      = 1
rho_rate = 0.9
mps_tol  = 1e-4
iter     = 0
eclipse  = False # Commented out in eom.py - too slow, need use a C++ eom and wrap in python)
while rho > 1e-4:

    if rho < 1e-3:
        mps_tol  = 1e-5
        rho_rate = 0.95
    # Method of Particular Solutions
    [soln,p0,mps_err] = mps_mee_ocp(tspan,p0,mee0,meef,rho,eclipse,mps_tol)
    print("rho:",rho,"mps err:",mps_err)
    rho  = rho*rho_rate

    # Update
    iter = iter + 1

# Final Solution
rho = rho/rho_rate # rho for converged solution
mee0_new = np.zeros(14)
mee0_new[0:7] = mee0[0:7]
mee0_new[7:14] = p0[0:7]

# Propagation Solution
sol   = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,mee0_new,method='LSODA',rtol=1e-13)
data  = np.array(soln.y,dtype=float).T
time  = np.array(soln.t,dtype=float).T

# Plot Solution
plot_mee_minfuel(time,data,rho,eclipse)

# Save Solution
df = pd.DataFrame(np.array([time,data]).T,columns=['time','trajectory'])
df.to_csv("output_soln.csv", index=False)
# Save Solution Parameters
df = pd.DataFrame(np.array([rho,eclipse,sc.Isp,sc.A,sc.eff]).T,columns=['rho','eclipse','Isp','Solar Array Area','Solar Array Efficiency'])
df.to_csv("output_params.csv", index=False)

# df = pd.read_csv("example.csv")
# print(df.values[:,0])
