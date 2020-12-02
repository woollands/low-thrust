import numpy as np
import const as cn
from scipy.integrate import solve_ivp as integrator
from eom import eom_mee_twobodyJ2_minfuel
from numba import njit
import sys

# @njit
def mee_ocp_central_difference(p0):
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
    states0 = np.array([1.8225634499541166, 0.725, 0.0, 0.06116262015048431, 0.0, 0, 100])
    statesf = np.array([6.610864583184714, 0.0, 0.0, 0.0, 0.0, 0.0])

    # Initialization
    IC        = np.zeros(14)
    Jac       = np.zeros((14,14))
    Del_plus  = np.zeros((7,7))
    Del_minus = np.zeros((7,7))
    DEL       = 1e-5*(p0/np.linalg.norm(p0))
    # DEL       = np.array([1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5])

    for tcnt in range(0,14):
        print(tcnt)
        if (tcnt < 7):
            IC[0:7]    = states0[0:7]
            IC[7:14]   = p0[0:7]
            IC[tcnt+7] = IC[tcnt+7]+DEL[tcnt]
        elif (tcnt >= 7):
            IC[0:7]  = states0[0:7]
            IC[7:14] = p0[0:7]
            IC[tcnt] = IC[tcnt]-DEL[tcnt-7]

        # Integrate Dynamics
        sol = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC,method='LSODA',rtol=1e-12)
        soln = sol.y[0:7,-1]
        if (tcnt < 7):
            Del_plus[0:7,tcnt] = soln[0:7]  # Free final mass & free final true longitude (angle)
        elif (tcnt >= 7):
            Del_minus[0:7,tcnt-7] = soln[0:7]  # Free final mass & free final true longitude (angle)

    # Compute Jacobian (central differences)
    for i in range(7):
        for j in range(7):
            Jac[i+7,j] = (Del_plus[i,j] - Del_minus[i,j])/(2*DEL[j])

    # print(Jac)
    # sys.exit()
    # MatInv = np.linalg.pinv(Del)
    # res    = statesf[0:5] - ref[0:5]     # Free final mass & free final true longitude (angle)
    # del_p0 = np.matmul(MatInv,res)

    return Jac
