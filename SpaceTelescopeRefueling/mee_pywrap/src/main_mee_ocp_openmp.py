import ctypes
import numpy as np
from functions import classical2mee
import sys
# Time
t0 = 0;
tf = 7121;
tspan = [t0,tf]

# Initial State
sma = 8000#24505
ecc = 0.1#0.725
inc = 10#7
Om  = 0
w   = 0
M   = 0
# m0  = 100
[p,f,g,h,k,L] = classical2mee(sma,ecc,np.deg2rad(inc),Om,w,M)
mee0 = [p, f, g, h, k, L]


# Final State
# sma = 42165
# ecc = 0
# inc = 0
# Om  = 0
# w   = 0
# M   = 0
# [p,f,g,h,k,L] = classical2mee(sma/6378.137,ecc,np.deg2rad(inc),Om,w,M)
# meef = [p, f, g, h, k, L]
# xT[0] = 7.055503292568207;
# xT[1] = 0.0;
# xT[2] = 0.0;
# xT[3] = 0.0;
# xT[4] = 0.0;
# xT[5] = 53.407075111026487;

# Set Initial Segment Scheme
# n = 2*np.pi/7;
# seg = n;
# TVEC_in = 0;

# u_inert = np.zeros((ind,3))

# sys.exit()
# for i = 1:length(Xp(:,6))
#     if Xp(i,6) > n
#         TVEC_in = [TVEC_in; Tp(i-1)];
#         n = n+seg;
#     end
# end
# TVEC_orig = TVEC_in;
# TVECprev_in = TVEC_in;
lib_mee_prop = ctypes.CDLL('./mee_prop.so')
lib_mee_prop.mee_prop.restype = None

lib_mee_prop.mee_prop.argtypes = (ctypes.c_int,
                                   ctypes.POINTER(ctypes.c_float),
                                   ctypes.c_int,
                                   ctypes.POINTER(ctypes.c_float))#,
                                   # ctypes.c_int,
                                   # ctypes.POINTER(ctypes.c_float))

def mee_prop(tspan,mee0):
    """Propagate MEEs using Picard-Chebyshev numerical integration."""
    meef  = np.zeros(6)
    sz1 = 2
    sz2 = 6
    sz3 = 6
    tspan = (ctypes.c_float * sz1)(*tspan)
    mee0  = (ctypes.c_float * sz2)(*mee0)
    meef  = (ctypes.c_float * sz3)(*meef)
    lib_mee_prop.mee_prop(sz1, tspan, sz2, mee0)#, sz3, meef)
    return meef

if __name__ == '__main__':
    statesf = mee_prop(tspan,mee0)
    # print('meef:\n', statesf[:])
    # print('\n')
