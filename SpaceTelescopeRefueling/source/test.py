import numpy as np
from scipy.optimize import minimize
from scipy.optimize import minimize
from scipy.integrate import solve_ivp as integrator
from eom import eom_mee_twobodyJ2_minfuel
import sys

import numpy as np
from scipy.optimize import minimize
# from numdifftools import Jacobian, Hessian

def fun(x, a):
    res = (x[0] - 1) **2 + (x[1] - a) **2
    # print(res)
    return res

def fun_der(x, a):
    # dx = 2 * (x[0] - 1)
    # dy = 2 * (x[1] - a)
    x_temp = np.zeros(2)
    d = 1e-5
    for i in range(2):
        x_temp[0] = x[0]
        x_temp[1] = x[1]
        if (i == 0):
            x_temp[i] = x[i]+d
            delp = fun(x_temp, a)
            x_temp[i] = x[i]-d
            delm = fun(x_temp, a)
            dx = (delp-delm)/2/d
        if (i == 1):
            x_temp[i] = x[i]+d
            delp = fun(x_temp, a)
            x_temp[i] = x[i]-d
            delm = fun(x_temp, a)
            dy = (delp-delm)/2/d

    return np.array([dx, dy])

def fun_hess(x, a):
    # dx = 2
    # dy = 2
    x_temp = np.zeros(2)
    d = 1e-5
    for i in range(2):
        x_temp[0] = x[0]
        x_temp[1] = x[1]
        if (i == 0):
            x_temp[i] = x[i]+d
            delp = fun_der(x_temp, a)
            x_temp[i] = x[i]-d
            delm = fun_der(x_temp, a)
            dx = (delp-delm)/2/d
        if (i == 1):
            x_temp[i] = x[i]+d
            delp = fun_der(x_temp, a)
            x_temp[i] = x[i]-d
            delm = fun_der(x_temp, a)
            dy = (delp-delm)/2/d

    M = np.zeros((2,2))
    M[0,:] = dx
    M[1,:] = dy
    return M

x0 = np.array([-100, 100]) # initial guess
a = 2.5

# res = minimize(fun, x0, args=(a,), method='dogleg', jac=fun_der, hess=fun_hess)
# print(res)

# sys.exit()
################################################################################
### The code below is an example problem copied from the scipy.optimize website
################################################################################

def rosen(x):
    """The Rosenbrock function"""
    res = sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)
    print(res)
    return res


def rosen_der(x):
    xm = x[1:-1]
    xm_m1 = x[:-2]
    xm_p1 = x[2:]
    der = np.zeros_like(x)
    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

def rosen_hess(x):
    x = np.asarray(x)
    H = np.diag(-400*x[:-1],1) - np.diag(400*x[:-1],-1)
    diagonal = np.zeros_like(x)
    diagonal[0] = 1200*x[0]**2-400*x[1]+2
    diagonal[-1] = 200
    diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]
    H = H + np.diag(diagonal)
    return H


x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
# print(x0)
# soln = minimize(rosen, x0, method='nelder-mead',
#                options={'xatol': 1e-8, 'disp': True})

# soln = minimize(rosen, x0, method='BFGS', jac=rosen_der,
#                options={'disp': True})

# soln = minimize(rosen, x0, method='Newton-CG',
#                jac=rosen_der, hess=rosen_hess,
#                options={'xtol': 1e-8, 'disp': True})

# soln = minimize(rosen, x0, method='dogleg',
#                jac=rosen_der, hess=rosen_hess,
#                options={'xtol': 1e-8, 'disp': True})

################################################################################
### The above code is an example problem copied from the scipy.optimize website
################################################################################

# def jacobian(p0):
#     # "Jacobian" of the residual function (scalar) with is required by some of the scipy optimization algorithms
#     # TEMP: For now I have copied these parameters into the function to make it more inline with the Rosenbrock
#     # function above, but eventually we'll need to pass them as input arguments
#     tspan = np.array([0,642.54])
#     rho = 1
#     eclipse = False
#     mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
#     meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])
#
#     # Initialization
#     IC1  = np.zeros(14)
#     IC2  = np.zeros(14)
#     grad = np.zeros(7)
#     DEL  = 1e-5*(p0/np.linalg.norm(p0))
#
#     for i in range(0,7):
#
#         IC1[0:7] = mee0
#         IC1[7:14] = p0
#         IC2[0:7] = mee0
#         IC2[7:14] = p0
#         IC1[i+7] = IC1[i+7]+DEL[i]
#         IC2[i+7] = IC2[i+7]-DEL[i]
#
#         # Integrate Dynamics
#         sol1 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC1,method='LSODA',rtol=1e-13)
#         sol2 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC2,method='LSODA',rtol=1e-13)
#
#         res1 = np.linalg.norm(meef[0:5] - sol1.y[0:5,-1])
#         res2 = np.linalg.norm(meef[0:5] - sol2.y[0:5,-1])
#
#         grad[i] = (res1-res2)/(2*DEL[i])
#
#     return grad


def fun_der_mee_ocp(x):

    # TEMP
    tspan = np.array([0,642.54])
    rho = 1
    eclipse = False
    mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
    meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])

    x_temp = np.zeros(7)
    d = 1e-5
    for i in range(7):
        x_temp[0] = x[0]
        x_temp[1] = x[1]
        x_temp[2] = x[2]
        x_temp[3] = x[3]
        x_temp[4] = x[4]
        x_temp[5] = x[5]
        x_temp[6] = x[6]
        if (i == 0):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp0 = (delp-delm)/2/d
        if (i == 1):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp1 = (delp-delm)/2/d
        if (i == 2):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp2 = (delp-delm)/2/d
        if (i == 3):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp3 = (delp-delm)/2/d
        if (i == 4):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp4 = (delp-delm)/2/d
        if (i == 5):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp5 = (delp-delm)/2/d
        if (i == 6):
            x_temp[i] = x[i]+d
            delp = residuals(x_temp)
            x_temp[i] = x[i]-d
            delm = residuals(x_temp)
            dp6 = (delp-delm)/2/d
    return np.array([dp0, dp1, dp2, dp3, dp4, dp5, dp6])

# def hessian(p0):
#     # "Hessian" of the residual function (scalar) with is required by some of the scipy optimization algorithms
#     # TEMP: For now I have copied these parameters into the function to make it more inline with the Rosenbrock
#     # function above, but eventually we'll need to pass them as input arguments
#     tspan = np.array([0,642.54])
#     rho = 1
#     eclipse = False
#     mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
#     meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])
#
#     # Initialization
#     IC1  = np.zeros(14)
#     IC2  = np.zeros(14)
#     IC3  = np.zeros(14)
#     IC4  = np.zeros(14)
#     grad = np.zeros(7)
#     hess = np.zeros((7,7))
#     DEL  = 1e-5*(p0/np.linalg.norm(p0))
#
#     for j in range(0,7):
#
#         IC1[0:7] = mee0
#         IC1[7:14] = p0
#         IC2[0:7] = mee0
#         IC2[7:14] = p0
#         IC1[j+7] = IC1[j+7]+DEL[j]
#         IC2[j+7] = IC2[j+7]-DEL[j]
#
#         for i in range(0,7):
#
#             IC3[0:7] = IC1[0:7]
#             IC3[7:14] = IC1[7:14]
#             IC4[0:7] = IC2[0:7]
#             IC4[7:14] = IC2[7:14]
#             IC3[i+7] = IC3[i+7]+DEL[i]
#             IC4[i+7] = IC4[i+7]-DEL[i]
#
#             # Integrate Dynamics
#             sol1 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC3,method='LSODA',rtol=1e-13)
#             sol2 = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC4,method='LSODA',rtol=1e-13)
#
#             res1 = np.linalg.norm(meef[0:5] - sol1.y[0:5,-1])
#             res2 = np.linalg.norm(meef[0:5] - sol2.y[0:5,-1])
#
#             hess[j,i] = (res1-res2)/(2*DEL[i])
#
#     return hess

def fun_hess_mee_ocp(x):

    # TEMP
    tspan = np.array([0,642.54])
    rho = 1
    eclipse = False
    mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
    meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])

    x_temp = np.zeros(7)
    d = 1e-5
    for i in range(7):
        x_temp[0] = x[0]
        x_temp[1] = x[1]
        x_temp[2] = x[2]
        x_temp[3] = x[3]
        x_temp[4] = x[4]
        x_temp[5] = x[5]
        x_temp[6] = x[6]
        if (i == 0):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp0 = (delp-delm)/2/d
        if (i == 1):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp1 = (delp-delm)/2/d
        if (i == 2):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp2 = (delp-delm)/2/d
        if (i == 3):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp3 = (delp-delm)/2/d
        if (i == 4):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp4 = (delp-delm)/2/d
        if (i == 5):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp5 = (delp-delm)/2/d
        if (i == 6):
            x_temp[i] = x[i]+d
            delp = fun_der_mee_ocp(x_temp)
            x_temp[i] = x[i]-d
            delm = fun_der_mee_ocp(x_temp)
            dp6 = (delp-delm)/2/d

    M = np.zeros((7,7))
    M[:,0] = dp0
    M[:,1] = dp1
    M[:,2] = dp2
    M[:,3] = dp3
    M[:,4] = dp4
    M[:,5] = dp5
    M[:,6] = dp6
    return M

def residuals(p0):
    # Residual function
    # TEMP: For now I have copied these parameters into the function to make it more inline with the Rosenbrock
    # function above, but eventually we'll need to pass them as input arguments
    tspan = np.array([0,642.54])
    rho = 1
    eclipse = False
    mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
    meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])

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


tspan = np.array([0,642.54])
rho = 1
eclipse = False
mee0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
meef = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])
p0 = np.zeros(7)
# p0 = np.random.rand(7)
p0[0] = -4.620949386264961
p0[1] = -12.037907266888872
p0[2] = 0.208961742408408
p0[3] = 5.946490544006020
p0[4] = 0.042440374825155
p0[5] = 0.002470486378745
p0[6] = 0.117893582793439

# res = minimize(residuals, p0, method='dogleg', jac=fun_der_mee_ocp, hess=fun_hess_mee_ocp)
# res = minimize(fun, x0, args=(a,), method='dogleg', jac=fun_der, hess=fun_hess)
# optres = minimize(residuals,p0,method='dogleg',jac=jacobian,shess=hessian,options={'xtol': 1e-8, 'disp': True})
# optres = minimize(residuals,p0,method='dogleg',jac=jacobian,hess=hessian,options={'xtol': 1e-8, 'disp': True})
optres = minimize(residuals, p0[0:7], jac='2-point', method = 'trust-constr')
print(optres)
# optres = minimize(residuals_mee_ocp,p0,method='Newton-CG',jac=mee_ocp_central_difference,options={'xtol': 1e-8, 'disp': True})
