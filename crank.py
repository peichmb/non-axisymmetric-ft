import numpy as np
from scipy.linalg.lapack import zgtsv
from params import *
import grid


dtype = 'complex'
tdSolve = zgtsv


def buildDiagonals(m, a0, a1):

    d = np.empty(nTheta, dtype=dtype)
    ud = np.empty(nTheta-1, dtype=dtype)
    ld = np.empty(nTheta-1, dtype=dtype)

    d[1:nTheta-1] = 1. + a0[1:nTheta-1] + a1[1:nTheta-1]
    ud[:] = -a1[:nTheta-1]
    ld[:] = -a0[1:]

    # Boundary conditions
    if m == 0:
        d[0] = 1. + a1[0]
        d[nTheta-1] = 1. + a0[nTheta-1]
    else:
        d[0] = 1. + 2.*a0[0] + a1[0]
        d[nTheta-1] = 1. + a0[nTheta-1] + 2.*a1[nTheta-1]

    return ld, d, ud


def buildRHS(m, a0, a1, b):

    rhs = np.empty(nTheta, dtype=dtype)
    rhs[1:nTheta-1] = a0[1:nTheta-1]*b[0:nTheta-2] \
                + (1. - a0[1:nTheta-1] - a1[1:nTheta-1])*b[1:nTheta-1] \
                + a1[1:nTheta-1]*b[2:nTheta]

    # Boundary conditions
    if m == 0:
        rhs[0] = (1. - a1[0])*b[0] + a1[0]*b[1]
        rhs[nTheta-1] = a0[nTheta-1]*b[nTheta-2] + (1.-a0[nTheta-1])*b[nTheta-1]
    else:
        rhs[0] = (1. - 2.*a0[0] - a1[0])*b[0] + a1[0]*b[1]
        rhs[nTheta-1] = a0[nTheta-1]*b[nTheta-2] + (1.-a0[nTheta-1]-2.*a1[nTheta-1])*b[nTheta-1]

    return rhs


def crankStep(b0, eta, dtime, m):

    a0 = dtime/2./grid.dTheta**2*grid.sinTheta_12[:-1]/grid.sinTheta*eta/rSun**2
    a1 = dtime/2./grid.dTheta**2*grid.sinTheta_12[1:]/grid.sinTheta*eta/rSun**2
    ld, d, ud = buildDiagonals(m, a0, a1)
    rhs = buildRHS(m, a0, a1, b0)

    du2, d, du, b1, info = tdSolve(ld, d, ud, rhs)

    return b1
