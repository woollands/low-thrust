import numpy as np
import const as cn
from eom import eom_mee_twobody_minfuel
from numba import njit

norm = np.linalg.norm

def mee_ocp_forward_difference():




    return


# @njit
def mee_ocp_central_difference(t,x,rho):

    delx = 1e-3*x/norm(x)

    dxdt_m = eom_mee_twobody_minfuel(t,x_m,rho)
    dxdt_p = eom_mee_twobody_minfuel(t,x_p,rho)

    for (j in range size(x)):
        jac[i,j] = (dxdt_p[i,j] - dxdt_m[i,j])/2/delx[j]


    return
