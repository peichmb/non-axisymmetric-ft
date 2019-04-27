import numpy as np
from params import *
import grid


def deriv_theta(f):
    
    f_theta = np.empty(nTheta, dtype='complex')
    # Centered finite differences
    f_theta[1:nTheta-1] = (f[2:nTheta] - f[0:nTheta-2])/2./grid.dTheta

    # Boundary conditions: f=0 @ poles
    f_theta[0] = (f[0] + f[1])/2./grid.dTheta
    f_theta[nTheta-1] = -(f[nTheta-1] + f[nTheta-2])/2./grid.dTheta

    return f_theta

#def deriv_theta(f):
#    
#    f_theta = np.empty(nTheta, dtype='complex')
#
#    # Centered finite differences
#    f_theta[2:nTheta-2] = (-f[4:nTheta]+8.*f[3:nTheta-1]-8.*f[1:nTheta-3]+f[nTheta-4])/12./grid.dTheta
#
#    # Boundary conditions: f=0 @ poles
#    f_theta[0] = (f[3]-6.*f[2]+18.*f[1]-7.*f[0])/12./grid.dTheta
#    f_theta[1] = (f[4]-6.*f[3]+18.*f[2]-10.*f[1]-3.*f[0])/12./grid.dTheta
#
#    f_theta[nTheta-2] = (3.*f[nTheta-1]+10.*f[nTheta-2]-18.*f[nTheta-3]+6.*f[nTheta-4]-f[nTheta-5])/12./grid.dTheta
#    f_theta[nTheta-1] = (7.*f[nTheta-1]-18.*f[nTheta-2]+6.*f[nTheta-3]-f[nTheta-4])/12./grid.dTheta
#
#    return f_theta

