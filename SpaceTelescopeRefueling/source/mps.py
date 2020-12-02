## @package mps
#  \brief     Method of particular solutions
#  \details   Solves two-point boundary value problems
#  \author    Robyn Woollands
#  \pre       numpy, const.py, spipy.integrate
#  \bug       No bugs known
#  \warning   May fail to converge for poor initial guess

import numpy as np
from scipy.integrate import solve_ivp as integrator
import const as cn
from eom import eom_twobody
from eom import eom_mee_twobodyJ2_minfuel
from numba import njit
import sys

__all__ = ["mps_twobody","mps_mee_ocp"]

################################################################################
def mps_twobody(tspan,states0,statesf):
    """
    Shooting method to solve Lambert's problem
    Parameters:
    ===========
    tspan   -- vector containing initial and final time
    states0 -- initial state vector (position & velocity)
    statesf -- final state vector (position & velocity)
    Returns:
    ========
    v0      -- initial velocity for transfer trajectory
    mps_err -- shooting method convergence error
    External:
    =========
    numpy, const.py, spipy.integrate
    """
    mps_err = 1
    mps_tol = 1e-13
    cnt     = 0
    IC      = np.zeros(6)
    Del     = np.zeros((3,3))
    DEL     = 1e-5
    v0      = states0[3:6]
    MatInv  = np.zeros((3,3))
    res     = np.zeros(3)
    del_v0  = np.zeros(3)
    while mps_err > mps_tol and cnt < 15:
        for tcnt in range(0,4):
            if tcnt == 0:
                IC[0:3] = states0[0:3]
                IC[3:6] = v0[0:3]
            elif tcnt > 0:
                IC[0:3] = states0[0:3]
                IC[3:6] = v0[0:3]
                IC[tcnt+2] = IC[tcnt+2]+DEL

            # Integrate Dynamics
            sol = integrator(eom_twobody,(tspan[0],tspan[-1]),IC,method='LSODA',rtol=1e-12)
            soln = sol.y[0:3,-1]
            if tcnt == 0:
                ref = soln
            elif tcnt > 0:
                Del[0:3,tcnt-1] = (soln[0:3]-ref[0:3])/DEL

        # Compute Updated Variables
        MatInv = np.linalg.inv(Del)
        res    = statesf[0:3] - ref[0:3]
        del_v0 = np.matmul(MatInv,res)
        v0     = v0 + del_v0

        # Compute Error
        mps_err = np.linalg.norm(res)/cn.DU
        # print(mps_err)

        # Counter
        cnt = cnt + 1

    return v0, mps_err
################################################################################
# @njit
def mps_mee_ocp(tspan,p0,states0,statesf,rho,eclipse,mps_tol):
    """
    Shooting method to solve min-fuel optimal control prolem
    Parameters:
    ===========
    tspan   -- vector containing initial and final time
    p0      -- initial costates guess
    states0 -- initial state vector (modified equinoctial elements)
    statesf -- final state vector (modified equinoctial elements)
    rho     -- switch smoothing parameter
    eclipse -- boolean (true or false)
    mps_tol -- user specified convergence tolerance
    Returns:
    ========
    traj    -- optimal trajectory
    p0      -- initial costate for optimal trajectory
    mps_err -- shooting method convergence error
    External:
    =========
    numpy, const.py, spipy.integrate
    """
    mps_err = 1
    cnt     = 0
    IC      = np.zeros(14)
    # Del     = np.zeros((5,7))
    # Del_plus  = np.zeros((5,7))
    # Del_minus = np.zeros((5,7))
    # Del       = np.zeros((7,7))
    Del       = np.zeros((6,7))
    Del_plus  = np.zeros((7,7))
    Del_minus = np.zeros((7,7))
    DEL     = np.array([1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5])
    # p0      = states0[7:14]
    MatInv  = np.zeros((7,6))
    # MatInv  = np.zeros((7,7))
    res     = np.zeros(6)
    del_p0  = np.zeros(7)
    while mps_err > mps_tol and cnt < 30:
        for tcnt in range(8): #range(0,15):
            if tcnt == 0:
                IC[0:7] = states0[0:7]
                IC[7:14] = p0[0:7]
            elif tcnt > 0 and tcnt < 8:
                IC[0:7] = states0[0:7]
                IC[7:14] = p0[0:7]
                IC[tcnt+6] = IC[tcnt+6]+DEL[tcnt-1]
            # elif tcnt > 7 and tcnt < 15:
            #     IC[0:7] = states0[0:7]
            #     IC[7:14] = p0[0:7]
            #     IC[tcnt+6-7] = IC[tcnt+6-7]-DEL[tcnt-1-7]

            # Integrate Dynamics
            sol = integrator(lambda t,y: eom_mee_twobodyJ2_minfuel(t,y,rho,eclipse),tspan,IC,method='LSODA',rtol=1e-13)
            # sol = integrator(eom_mee_twobody_minfuel,(tspan[0],tspan[-1]),IC,method='LSODA',rtol=1e-12)
            soln = sol.y[:,-1]
            if tcnt == 0:
                ref  = soln
                traj = sol
            elif tcnt > 0 and tcnt < 8:
                Del[0:5,tcnt-1] = (soln[0:5]-ref[0:5])/1e-5  # Free final mass & free final true longitude (angle)
                Del[5,tcnt-1] = (soln[13]-ref[13])/1e-5  # Free final mass & free final true longitude (angle)
            #     Del_plus[0:5,tcnt-1] = soln[0:5]  # Free final mass & free final true longitude (angle)
            #     Del_plus[5:7,tcnt-1] = soln[12:14]  # Free final mass & free final true longitude (angle)
            # elif tcnt > 7 and tcnt < 15:
            #     Del_minus[0:5,tcnt-1-7] = soln[0:5]  # Free final mass & free final true longitude (angle)
            #     Del_minus[5:7,tcnt-1-7] = soln[12:14]  # Free final mass & free final true longitude (angle)

        # Compute Updated Variables
        # for i in range(7): #range(5):
        #     for j in range(7):
        #         # Del[i,j] = (Del_plus[i,j] - Del_minus[i,j])/(2*DEL[j])
        #         Del[i,j] = Del[i,j]/DEL[j]

        # print(Del)
        MatInv   = np.linalg.pinv(Del)
        res[0:5] = statesf[0:5] - ref[0:5]     # Free final mass & free final true longitude (angle)
        res[5] = 0 - ref[13]     # Free final mass & free final true longitude (angle)
        del_p0   = np.matmul(MatInv,res)
        p0       = p0 + del_p0
        # DEL      = 1e-5*(p0/np.linalg.norm(p0))

        # Compute Error
        mps_err = np.linalg.norm(res)
        print(mps_err)

        # Counter
        cnt = cnt + 1

    return traj, p0, mps_err
################################################################################
