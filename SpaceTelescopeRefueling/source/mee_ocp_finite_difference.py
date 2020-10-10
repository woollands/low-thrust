import numpy as np
import const as cn
from eom import eom_mee_twobodyJ2_minfuel
from numba import njit

# @njit
def mee_ocp_central_difference(states0,tspan,p0,rho,eclipse):
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

    # Sizes
    sz1 = np.size(states0)
    sz2 = sz1 - 2 # Number of fixed final states

    # Initialization
    IC        = np.zeros(2*sz1)
    Jac       = np.zeros((sz2,sz1))
    Del_plus  = np.zeros((sz2,sz1))
    Del_minus = np.zeros((sz2,sz1))
    DEL       = 1e-5*(p0/np.linalg.norm(p0))

    for tcnt in range(0,2*sz1):
        if (tcnt < sz1):
            IC[0:sz1]    = states0[0:sz1]
            IC[0:sz1]    = p0[0:sz1]
            IC[tcnt+sz1] = IC[tcnt+sz1]+DEL[tcnt]
        elif (tcnt >= sz1):
            IC[0:sz1]     = states0[0:sz1]
            IC[sz1:2*sz1] = p0[0:sz1]
            IC[tcnt]      = IC[tcnt]-DEL[tcnt-sz1]

        # Integrate Dynamics
        sol = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC,method='LSODA',rtol=1e-13)
        soln = sol.y[0:sz1,-1]
        if (tcnt < sz1):
            Del_plus[0:sz2,tcnt-1] = soln[0:sz2]  # Free final mass & free final true longitude (angle)
        elif (tcnt >= sz1):
            Del_minus[0:sz2,tcnt-sz1] = soln[0:sz2]  # Free final mass & free final true longitude (angle)

    # Compute Jacobian (central differences)
    for i in range(sz2):
        for j in range(sz1):
            Jac[i,j] = (Del_plus[i,j] - Del_minus[i,j])/(2*DEL[j])

    return Jac
