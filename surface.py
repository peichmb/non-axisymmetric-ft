# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- surface.py: integrates the equation describing the evolution of
# ---   the surface magnetic field (br).
# --- ----------------------------------------------------------------

import numpy as np
from params import *
import grid
import flows
import crank
import derivatives


def rhs_advTheta(br_m):

    rhs = -1./grid.sinTheta/rSun*derivatives.deriv_theta(flows.uTheta_rSun*br_m*grid.sinTheta)

    return rhs


def rhs_latShear(m, br_m):

    rhs = -flows.omega_rSun*1j*m*br_m

    return rhs


def rhs_phiDiff(m, br_m):

    rhs = eta_surf/rSun**2/grid.sin2theta*(-m**2)*br_m 

    return rhs


def advance(m, br_m, dtime):

    br_m_1 = br_m + rhs_advTheta(br_m)*dtime
    br_m_1 = br_m_1 + rhs_latShear(br_m_1, m)*dtime
    br_m_1 = crank.crankStep(br_m_1, eta_surf, dtime, m)
    br_m_1 = br_m_1*np.exp(-eta_surf/rSun**2/grid.sin2theta*m**2*dtime)

    return br_m_1


