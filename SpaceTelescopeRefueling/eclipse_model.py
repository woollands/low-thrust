import numpy as np
import const as cn
import spiceypy as sp
import math
from functions import mee2rv
import sys

# Shortcuts
asin = math.asin
acos = math.acos
cos  = np.cos
tanh = np.tanh
norm = np.linalg.norm
sqrt = np.sqrt
dot  = np.dot

# Eclipse model - Grazing Goat
def eclipse_grazing_goat(t,X,rho):
    # Convert time back to units of seconds
    t = t*cn.TU
    # Convert to Cartesian
    [r,v] = mee2rv(X[0],X[1],X[2],X[3],X[4],X[5])
    rSC = r*cn.DU
    # Sun's state w.r.t. Earth
    [states,lt] = sp.spkezr('SUN', t, 'ECLIPJ2000', 'LT+S', 'EARTH');
    rS    = states[0:3]     # Sun's position vector
    rS_SC = rS - rSC        # Spacecraft's position w.r.t. Sun
    rB_SC = -rSC            # Spacecraft's position w.r.t. Body (Earth)
    Xp    = cn.Req/(cn.Req+cn.RS)*rS
    alpha = asin(cn.Req/norm(Xp))
    ang   = acos(dot(-(rSC-Xp),rS)/(norm(rSC-Xp)*norm(rS)))
    TF    = 1
    PA    = cn.Power
    zetaE = 1
    # Eclipse
    if rho < 1:
        if (dot(-(rSC-Xp),rS) >= norm(rS-Xp)*cos(alpha)):
            Rapp = asin(cn.RS/norm(rS_SC))    # Sun's apparent radius
            rapp = asin(cn.Req/norm(rB_SC))   # Earth's apparent radius
            d = acos(dot(rB_SC,rS_SC.T)/(norm(rB_SC)*norm(rS_SC)))   # apparent distance between disk centers
            # Hyperbolic Tangent Smoothing
            gammaE = (d-rapp+Rapp)/(rapp+Rapp)  # no eclipse (g > 0), eclipse (g < 0)
            zetaE  = 0.5*(1+tanh(gammaE/rho)) # smoothing parameter
            d2 = d**2
            r2 = rapp**2
            R2 = Rapp**2
            # Grazing Goat Model
            if (d < 1e-10):
                if (rapp >= Rapp): # total eclipse
                    TF = 0
                    PA   = zetaE*cn.Power # Power Available
                else: # transit
                    A    = pi*(Rapp**2 - rapp**2)
                    Atot = pi*Rapp**2
                    TF   = 1 - A/Atot   # Transit Factor
                    PA   = zetaE*cn.Power # Power Available
            else:
                # var1 = (d2 +r2 - R2)/(2*d*rapp)
                # var2 = (d2 - r2 + R2)/(2*d*Rapp)
                # var3 = (d+rapp-Rapp)*(d-rapp+Rapp)*(-d+rapp+Rapp)*(d+rapp+Rapp)
                # if (var1 > 1):
                #     var1 = 1
                # elif (var1 < -1):
                #     var1 = -1
                # if (var2 > 1):
                #     var2 = 1
                # elif (var2 < -1):
                #     var2 = -1
                # if (var3 < 0):
                #     var3 = 0
                # A    = r2 * acos((var1)) + R2 * acos((var2)) - 0.5*sqrt((var3))
                # Atot = pi*Rapp**2
                # TF   = 1 - A/Atot  # Transit Factor
                PA   = zetaE*cn.Power # Power Available
    return TF, PA, zetaE

# Eclipse Model - Cylindrical
def eclipse_cylindrical(t,X,rho):
    [r,v] = mee2rv(X[0],X[1],X[2],X[3],X[4],X[5])
    gamma = norm(temp[1:2]) - 1   # Assume Sun is located along positive x-axis
    zeta  = 0.5*(1+tanh(gamma/rho))
    if (rSC[0] < 0):
        PA = cn.Power*zeta
        TF = zeta
    else:
        PA = cn.Power
        TF = 1
    return PA
