# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- flows.py: contains the definition of the large scale flows:
# ---   shear by differential rotation and meridional circulation.
# --- ----------------------------------------------------------------

import numpy as np
from params import *
import grid

# Effective theta-velocity for the advection term (toroidal field)
uTheta = mflow_return*np.sin(2.4*grid.lat)
uTheta[np.where(grid.lat>75.*pi/180.)] = 0.
uTheta[np.where(grid.lat<-75.*pi/180.)] = 0.
# Mflow goes up to the pole
uTheta = mflow_return*np.sin(2.*grid.lat)

# Surface meridional flow
uTheta_rSun = -mflow_surf*np.sin(2.4*grid.lat)
uTheta_rSun[np.where(grid.lat>75.*pi/180.)] = 0.
uTheta_rSun[np.where(grid.lat<-75.*pi/180.)] = 0.
# Mflow goes up to the pole
uTheta_rSun = -mflow_surf*np.sin(2.*grid.lat)

# Differential rotation
dd2rs = (pi/180.)/(24.*3600.) # degrees/day to rad/sec
if dRotModel == 1:

    omega_rSun = 14.306 - 1.98*grid.cos2theta - 2.14*grid.cos4theta
    omega_rNSSL = omega_rSun + 0.53
    dOmega_dTheta_rNSSL = 1.98*2.*grid.cosTheta*grid.sinTheta \
                         +2.14*4.*grid.cos3theta*grid.sinTheta
    dOmega_dTheta_rNSSL = dOmega_dTheta_rNSSL*dd2rs
    omega_rNSSL = (omega_rNSSL-360./(27.2753-1.77))*dd2rs # -synodic... ; rad/s
    omega_rSun = (omega_rSun-360./(27.2753-1.77))*dd2rs # -synodic... ; rad/s
    deltaOmega = omega_rSun - omega_rNSSL

if dRotModel == 2:

    omega_rSun = 14.437 - 1.48*grid.cos2theta - 2.99*grid.cos4theta #From Hathaway & Rightmire 2011
    dOmega_dTheta_rNSSL = 1.48*2.*grid.cosTheta*grid.sinTheta + 2.99*4.*grid.cos3theta*grid.sinTheta

    omega_rSun_helio = 14.1772 - 1.59252*grid.cos2theta - 2.61274*grid.cos4theta # From Schou et al 1998
    omega_rNSSL= omega_rSun_helio + 0.53
    #omega_rSun = 14.306 - 1.98*grid.cos2theta - 2.14*grid.cos4theta
    #omega_rNSSL = omega_rSun + 0.53
    #dOmega_dTheta_rNSSL = 1.98*2.*grid.cosTheta*grid.sinTheta \
    #                     +2.14*4.*grid.cos3theta*grid.sinTheta
    dOmega_dTheta_rNSSL = dOmega_dTheta_rNSSL*dd2rs
    omega_rNSSL = (omega_rNSSL-360./(27.2753-1.77))*dd2rs # -synodic... ; rad/s
    omega_rSun = (omega_rSun-360./(27.2753-1.77))*dd2rs # -synodic... ; rad/s
    #deltaOmega = omega_rSun - omega_rNSSL




#import matplotlib.pyplot as plt
#plt.plot(grid.cosTheta,dOmega_dTheta_rNSSL )

#plt.show()
#exit()

