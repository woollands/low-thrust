## @package plot_routines
#  \brief     Plot solution
#  \details   Plots the states, costates switch function, thrust profile and
#  3D trajectory
#  \author    Robyn Woollands
#  \pre       numpy, pandas, mpl_toolkits, matplotlib, PIL
#  \bug       No bugs known

import pandas as pd
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from mpl_toolkits.basemap import Basemap
from functions import mee2rv
from functions import thrust_angle

__all__ = ["plot_mee_minfuel","drawEarth"]

################################################################################
def plot_mee_minfuel(t,data,rho,eclipse):
    """
    Plot solution
    Parameters:
    ===========
    t       -- time
    data    -- state and costate time history
    rho     --
    eclipse --
    Returns:
    ========
    External:
    =========
    numpy, pandas, mpl_toolkits, matplotlib
    """
    p = data[:,0]
    f = data[:,1]
    g = data[:,2]
    h = data[:,3]
    k = data[:,4]
    L = data[:,5]
    m = data[:,6]
    plam = data[:,7]
    flam = data[:,8]
    glam = data[:,9]
    hlam = data[:,10]
    klam = data[:,11]
    Llam = data[:,12]
    mlam = data[:,13]

    [r,v,u_inert,u_lvlh,S,F,Pa,delta,zeta] = thrust_angle(data,rho,eclipse)

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('States 1')
    ax1.plot(t,p)
    ax1.set_ylabel('p (DU)')
    ax2.plot(t,f)
    ax2.set_ylabel('f')
    ax3.plot(t,g)
    ax3.set_xlabel('Time (TU)')
    ax3.set_ylabel('g')

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('States 2')
    ax1.plot(t,h)
    ax1.set_ylabel('h')
    ax2.plot(t,k)
    ax2.set_ylabel('k')
    ax3.plot(t,L)
    ax3.set_xlabel('Time (TU)')
    ax3.set_ylabel('L')

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Costates 1')
    ax1.plot(t,plam)
    ax1.set_ylabel(r'$\lambda_{\mathrm{p}}$')
    ax2.plot(t,flam)
    ax2.set_ylabel(r'$\lambda_{\mathrm{f}}$')
    ax3.plot(t,glam)
    ax3.set_xlabel('Time (TU)')
    ax3.set_ylabel(r'$\lambda_{\mathrm{g}}$')

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Costates 2')
    ax1.plot(t,hlam)
    ax1.set_ylabel(r'$\lambda_{\mathrm{h}}$')
    ax2.plot(t,klam)
    ax2.set_ylabel(r'$\lambda_{\mathrm{k}}$')
    ax3.plot(t,Llam)
    ax3.set_xlabel('Time (TU)')
    ax3.set_ylabel(r'$\lambda_{\mathrm{L}}$')

    fig, (ax1, ax2) = plt.subplots(2)
    fig.suptitle('Mass')
    ax1.plot(t,m)
    ax1.set_ylabel('Mass (kg)')
    ax2.plot(t,mlam)
    ax2.set_ylabel(r'$\lambda_{\mathrm{m}}$')
    ax2.set_xlabel('Time (TU)')
    ax1.grid()
    ax2.grid()

    # fig, (ax1, ax2) = plt.subplots(2)
    # fig.suptitle('Power')
    # ax1.plot(t,Pa)
    # ax1.set_ylabel('Available Power')
    # ax2.plot(t,delta)
    # ax2.plot(t,zeta)
    # ax2.plot(t,delta*zeta)
    # ax2.legend([r'$\delta$',r'$\zeta$',r'$\delta$ * $\zeta$'])
    # ax2.set_ylabel('Continuation Parameters')
    # ax2.set_xlabel('Time (TU)')
    # ax1.grid()
    # ax2.grid()

    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.suptitle('Power and Thrust')
    ax1.plot(t,S)
    ax1.set_ylabel('Switch Function')
    ax2.plot(t,Pa)
    ax2.set_ylabel('Power (W)')
    ax3.plot(t,F)
    ax3.set_ylabel('Thrust (N)')
    ax3.set_xlabel('Time (TU)')
    ax1.grid()
    ax2.grid()
    ax3.grid()

    [x, y, z, bm] = drawEarth(1)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(r[:,0],r[:,1],r[:,2])
    ax.quiver(r[:,0],r[:,1],r[:,2],u_lvlh[:,0],u_lvlh[:,1],u_lvlh[:,2],length=0.2,colors='red')
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, alpha=0.5, facecolors=bm)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_zlabel('z (km)')
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,6)
    ax.set_zlim(-6,6)
    # axis.set_clip_on(self,b)

    plt.show()

################################################################################
def drawEarth(Radius):
    """
    Computes meshgrid for plotting the Earth
    Parameters:
    ===========
    Radius -- radius of Earth
    Returns:
    ========
    x   -- x position
    y   -- y position
    z   -- z position
    img -- RGB values
    External:
    =========
    numpy, PIL, blue_marble.jpg
    """
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
################################################################################
