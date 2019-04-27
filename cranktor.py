import numpy as np
from scipy.linalg.lapack import zgtsv
from params import *
import grid


dtype = 'complex'
tdSolve = zgtsv


def buildDiagonals(m, a0, a1, dtime):

    d = np.empty(nTheta, dtype=dtype)
    ud = np.empty(nTheta-1, dtype=dtype)
    ld = np.empty(nTheta-1, dtype=dtype)

    d[1:nTheta-1] = 1. + a0[1:nTheta-1] + a1[1:nTheta-1]
    ud[:] = -a1[:nTheta-1]
    ld[:] = -a0[1:]

    # Boundary conditions
    d[0] = 1.-dtime*eta0/rSun**2/2.*(grid.cosTheta[0]/grid.sinTheta[0]/2./grid.dTheta \
                                     -3./grid.dTheta**2)

    d[nTheta-1] = 1.-dtime*eta0/rSun**2/2.*(-grid.cosTheta[nTheta-1]/grid.sinTheta[nTheta-1]/2./grid.dTheta \
                                            -3./grid.dTheta**2)

    ud[0] = dtime*eta0/rSun**2/2.*(grid.cosTheta[0]/grid.sinTheta[0]/2./grid.dTheta \
                                   + 1./grid.dTheta**2)

    ld[nTheta-2] = dtime*eta0/rSun**2/2.*(-grid.cosTheta[nTheta-1]/grid.sinTheta[nTheta-1]/2./grid.dTheta \
                                          + 1./grid.dTheta**2)

    return ld, d, ud


def buildRHS(m, a0, a1, b, dtime):

    rhs = np.empty(nTheta, dtype=dtype)
    rhs[1:nTheta-1] = a0[1:nTheta-1]*b[0:nTheta-2] \
                + (1. - a0[1:nTheta-1] - a1[1:nTheta-1])*b[1:nTheta-1] \
                + a1[1:nTheta-1]*b[2:nTheta]

    # Boundary conditions
    rhs[0] = (1.+dtime*eta0/rSun**2/2.*(grid.cosTheta[0]/grid.sinTheta[0]/2./grid.dTheta \
                                     -3./grid.dTheta**2))*b[0] \
           + dtime*eta0/rSun**2/2.*(grid.cosTheta[0]/grid.sinTheta[0]/2./grid.dTheta \
                                   + 1./grid.dTheta**2)*b[1]

    rhs[nTheta-1] = dtime*eta0/rSun**2/2.*(-grid.cosTheta[nTheta-1]/grid.sinTheta[nTheta-1]/2./grid.dTheta \
                                          + 1./grid.dTheta**2)*b[nTheta-2] \
                  + (1.-dtime*eta0/rSun**2/2.*(-grid.cosTheta[nTheta-1]/grid.sinTheta[nTheta-1]/2./grid.dTheta \
                                            -3./grid.dTheta**2))*b[nTheta-1]

    return rhs


def crankStep(b0, eta, dtime, m):

    a0 = dtime/2./grid.dTheta**2*grid.sinTheta_12[:-1]/grid.sinTheta*eta/rSun**2
    a1 = dtime/2./grid.dTheta**2*grid.sinTheta_12[1:]/grid.sinTheta*eta/rSun**2
    ld, d, ud = buildDiagonals(m, a0, a1, dtime)
    rhs = buildRHS(m, a0, a1, b0, dtime)

    du2, d, du, b1, info = tdSolve(ld, d, ud, rhs)

    return b1
