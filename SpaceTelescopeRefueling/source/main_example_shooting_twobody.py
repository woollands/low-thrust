## @package main_example_shooting_twobody
#  \brief     Main file for two-body orbit transfer problem setup
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
from matplotlib import rc
import const as cn
import sys
from mps import mps_twobody
from eom import eom_twobody
# import plotly.express as px

def residuals(v0,tspan,states0,statesf):
    """
    Computes error between current final state and desired final state for the
    Keplerian Lambert type problem. This is a numerical method not an analytic
    Lambert solver.
    Parameters:
    ===========
    v0      -- current initial velocity
    tspan   -- vector containing initial and final time
    states0 -- initial states
    statesf -- final states
    Returns:
    ========
    res -- residual
    External:
    =========
    numpy, scipy.integrate
    """
    states_in = np.zeros(6)
    states_in[0:3] = states0[0:3]
    states_in[3:6] = v0[0:3]
    sol = integrator(eom_twobody,(tspan[0],tspan[-1]),states_in,method='LSODA',rtol=1e-12)
    res = np.linalg.norm(statesf[0:3] - sol.y[0:3,-1])/cn.DU
    return res

# Initial & Final Conditions
sma     = 8000
ecc     = 0
rp      = sma*(1-ecc)
P       = 2*np.pi*np.sqrt(sma**3 / cn.mu)
tspan   = np.array([0, 0.8*P])
elm0    = np.array([rp, ecc, np.deg2rad(10), np.deg2rad(0), np.deg2rad(0), np.deg2rad(0), 0, cn.mu])
states0 = sp.conics(elm0,0)
sma     = sma + 100;
rp      = sma*(1-ecc)
elm1    = np.array([rp, ecc, np.deg2rad(10), np.deg2rad(0), np.deg2rad(0), np.deg2rad(0), 0, cn.mu])
states1 = sp.conics(elm1,0)
statesf = sp.prop2b(cn.mu,states1,tspan[-1])

# Method of Particular Solutions Shooting Method
[v0,mps_err] = mps_twobody(tspan,states0,statesf)
print("mps err:",mps_err)
states0_new = np.zeros(6)
states0_new[0:3] = states0[0:3]
states0_new[3:6] = v0[0:3]

# Scipy's Optimize
v0_guess = states0[3:6]
optres = minimize(residuals,v0_guess[0:3],args=(tspan,states0,statesf),method='Nelder-Mead',tol=1e-13)
states0_new = np.zeros(6)
states0_new[0:3] = states0[0:3]
states0_new[3:6] = optres.x[0:3]
print("opt_err:",optres.fun)

# Check Solution
sol = integrator(eom_twobody,(tspan[0],tspan[-1]),states0_new,method='LSODA',rtol=1e-12)

# Plot
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(sol.y[0,:],sol.y[1,:],sol.y[2,:])
ax.scatter3D(states0[0],states0[1],states0[2],'r')
ax.scatter3D(statesf[0],statesf[1],statesf[2],'g')
plt.show()

# import pandas as pd
# statesdf=pd.DataFrame([states0,statesf] ,columns=['x','y','z','vx','vy','vz'])
# statesdf['cat']=['s','f']
# soldf=pd.DataFrame(sol.y.T ,columns=['x','y','z','vx','vy','vz'])
# soldf['cat']='o'
# soldf2=pd.concat([soldf,statesdf]).reset_index(drop=True)
# print(soldf2)
# fig = px.scatter_3d(soldf2,x='x', y='y', z='z',color='cat')
# fig.show()
