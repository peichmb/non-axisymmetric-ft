# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- grid.py: This module contains the grid and associated parameters
# --- ----------------------------------------------------------------

import numpy as np
from params import *

dPhi = 2.*pi/nPhi
phi = np.linspace(dPhi/2., 2.*pi-dPhi/2., nPhi)

# Grid is equally spaced in theta
dTheta = pi/nTheta
theta = np.linspace(dTheta/2., pi-dTheta/2., nTheta)
lat = pi/2.-theta

sinTheta = np.sin(theta)
cosTheta = np.cos(theta)
sin2theta = sinTheta**2
sin3theta = sinTheta**3
sin4theta = sinTheta**4
cos2theta = cosTheta**2
cos3theta = cosTheta**3
cos4theta = cosTheta**4

# These refer to the cell boundaries.
# Needed for the Crank-Nicolson scheme.
theta_12 = np.linspace(0, pi, nTheta+1)
sinTheta_12 = np.sin(theta_12)
sin2Theta_12 = np.sin(theta_12)**2

# Azimuthal orders
m = np.fft.rfftfreq(nPhi)*nPhi
