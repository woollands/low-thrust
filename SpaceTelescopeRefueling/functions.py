import numpy as np
from numpy import linalg as LA
from PIL import Image
import spacecraft_params as sc
import const as cn
from numba import njit

# Shortcuts
cos    = np.cos
sin    = np.sin
tan    = np.tan
sqrt   = np.sqrt
tanh   = np.tanh
cross  = np.cross
matmul = np.matmul
norm   = LA.norm


__all__ = ["mee2rv", "inertial2radial", "thrust_angle", "drawEarth"]

################################################################################
# @njit
def mee2rv(p,f,g,h,k,L):
    """ Convert modified equinoctial elements to cartesian coordinates
    Robyn Woollands 04/13/2020
    Parameters:
    p  -- semilatus rectum
    f  -- e*cos(w+Omega)
    g  -- e*sin(w+Omega)
    h  -- tan(i/2)*cos(Omega)
    k  -- tan(i/2)*sin(Omega)
    L  -- true longitude
    mu -- standard gravitational parameter
    Returns:
    r -- position
    v -- velocity
    """
    # Common terms
    alpha2 = h**2 - k**2
    tani2s = h**2 + k**2
    s2     = 1 + tani2s
    cosL   = cos(L)
    sinL   = sin(L)
    hk2    = 2*h*k
    sq     = sqrt(cn.mu/p)
    # Radius
    radius = p/(1 + f*cosL + g*sinL)
    # Position and Velocity
    r = np.zeros((np.size(p),3))
    v = np.zeros((np.size(p),3))
    r[:,0] = radius*(cosL + alpha2*cosL + hk2*sinL)/s2
    r[:,1] = radius*(sinL - alpha2*sinL + hk2*cosL)/s2
    r[:,2] = 2*radius*(h*sinL - k*cosL)/s2
    v[:,0] = -sq*(sinL + alpha2*sinL - hk2*cosL + g - f*hk2 + alpha2*g)/s2
    v[:,1] = -sq*(-cosL + alpha2*cosL + hk2*sinL - f + g*hk2 + alpha2*f)/s2
    v[:,2] = 2*sq*(h*cosL + k*sinL + f*h + g*k)/s2
    return r, v
################################################################################
# @njit
def classical2mee(a,e,inc,Om,w,nu):
    p = a*(1 - e**2);
    f = e*cos(w + Om);
    g = e*sin(w + Om);
    h = tan(inc/2)*cos(Om);
    k = tan(inc/2)*sin(Om);
    L = Om + w + nu;
    return p, f, g, h, k, L
################################################################################
# @njit
def inertial2radial(r,v):
    # Radial in x- and z- direction
    hvec = cross(r,v);
    h    = norm(hvec);
    xrdl = r/norm(r);
    zrdl = hvec/h;
    # Radial in y-direction
    yrdl = cross(zrdl,xrdl);
    return np.array([xrdl,yrdl,zrdl])
################################################################################
# @njit
def switch_function(data,rho):
    # Assign Variables
    p       = data[:,0]
    f       = data[:,1]
    g       = data[:,2]
    h       = data[:,3]
    k       = data[:,4]
    L       = data[:,5]
    m       = data[:,6]
    plam    = data[:,7]
    flam    = data[:,8]
    glam    = data[:,9]
    hlam    = data[:,10]
    klam    = data[:,11]
    Llam    = data[:,12]
    mlam    = data[:,13]
    # Initialization
    ind     = np.size(data,0)
    S       = np.zeros(ind)
    delta   = np.zeros(ind)
    BTL     = []
    for i in range(ind):
        # Common terms
        SinL = sin(L[i])
        CosL = cos(L[i])
        q    = 1+f[i]*CosL+g[i]*SinL
        s    = 1+h[i]**2+k[i]**2
        C1   = sqrt(p[i])
        C2   = 1/q
        C3   = h[i]*SinL-k[i]*CosL
        B = np.array([[0,2*p[i]*C2*C1,0],
        [C1*SinL,C1*C2*((q+1)*CosL+f[i]),-C1*(g[i]/q)*C3],
        [-C1*CosL,C1*C2*((q+1)*SinL+g[i]),C1*(f[i]/q)*C3],
        [0,0,C1*s*CosL*C2/2],
        [0,0,C1*s*SinL*C2/2],
        [0,0,C1*C2*C3],
        [0,0,0]])
        lam = np.array([plam[i],flam[i],glam[i],hlam[i],klam[i],Llam[i],mlam[i]])
        BTL.append(matmul(B.T,lam))
        S[i]     = sc.c*sc.si2can*norm(BTL[i])/m[i]+mlam[i]-1
        delta[i] = 0.5*(1+tanh(S[i]/rho))
    return S, delta, BTL
################################################################################
# @njit
def thrust_angle(data,rho,eclipse):
    # Assign Variables
    p       = data[:,0]
    f       = data[:,1]
    g       = data[:,2]
    h       = data[:,3]
    k       = data[:,4]
    L       = data[:,5]
    m       = data[:,6]
    plam    = data[:,7]
    flam    = data[:,8]
    glam    = data[:,9]
    hlam    = data[:,10]
    klam    = data[:,11]
    Llam    = data[:,12]
    mlam    = data[:,13]
    # Initialization
    ind     = np.size(data,0)
    F       = np.zeros(ind)
    zeta    = np.zeros(ind)
    Pa      = np.zeros(ind)
    u_inert = np.zeros((ind,3))
    u_lvlh  = np.zeros((ind,3))
    # Switch Function
    [S,delta,BTL] = switch_function(data,rho)
    [r,v] = mee2rv(p,f,g,h,k,L)
    for i in range(ind):
        # Eclipse model
        if (eclipse):
            if (r[0] < 0): # Assume Sun is located along positive x-axis
                gamma   = norm(r[1:2]) - cn.Req/cn.DU
                zeta[i] = 0.5*(1.0+tanh(gamma/rho))
                Pa[i]   = zeta[i]*cn.Power
                Thr     = Pa[i]*sc.A*sc.eta/sc.Isp/cn.g0
                F[i]    = Thr*delta[i]
                u_inert[i,:] = -BTL[i]/norm(BTL[i])*delta[i]*zeta[i]
            else:
                Pa[i] = cn.Power
                Thr   = Pa[i]*sc.A*sc.eta/sc.Isp/cn.g0
                F[i]  = Thr*delta[i]
                u_inert[i,:] = -BTL[i]/norm(BTL[i])*delta[i]
        else:
            Pa[i] = cn.Power
            Thr   = cn.P*sc.A*sc.eta/sc.Isp/cn.g0
            F[i]  = Thr*delta[i]
            u_inert[i,:] = -BTL[i]/norm(BTL[i])*delta[i]

        M = inertial2radial(r[i,:],v[i,:])
        u_lvlh[i,:] = matmul(M,u_inert[i,:].T).T
    return r, v, u_inert, u_lvlh, S, F, Pa, delta, zeta
################################################################################
def drawEarth(Radius):

    # Create a sphere with earths surface texture

    # Load texture
    #response = requests.get('http://www.johnstonsarchive.net/spaceart/cmaps/earthmap.jpg')

    #img = Image.open(StringIO(response.content))
    img = Image.open('blue_marble.jpg')

    # Rescale RGB values
    img = np.array(img.resize([int(d/4) for d in img.size]))/256

    # Image coordinates
    lons = np.linspace(-180, 180, img.shape[1]) * np.pi/180
    lats = np.linspace(-90, 90, img.shape[0])[::-1] * np.pi/180

    x = Radius*np.outer(np.cos(lons), np.cos(lats)).T
    y = Radius*np.outer(np.sin(lons), np.cos(lats)).T
    z = Radius*np.outer(np.ones(np.size(lons)), np.sin(lats)).T

    # Alternatively, create a simple sphere object (faster)
    # pi = np.pi
    # phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
    # x = Radius*np.sin(phi)*np.cos(theta)
    # y = Radius*np.sin(phi)*np.sin(theta)
    # z = Radius*np.cos(phi)

    return x, y, z, img
