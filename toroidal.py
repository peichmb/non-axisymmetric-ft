# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- toroidal.py: integrates the equation describing the evolution of
# ---   the radially ingegrated toroidal flux (bTor).
# --- ----------------------------------------------------------------

from params import *
import grid
import flows
import cranktor
import derivatives

mMax = 32 

def calcThetaComp(m, br_m, bTor_m):
    
    # Update this routine to use the trapezoidal rule
    f1 = grid.sinTheta*br_m*grid.dTheta
    f2 = bTor_m*grid.dTheta
    bTht_m = - rSun**2*(np.cumsum(f1)-0.5*f1) \
             - 1j*m*(np.cumsum(f2)-0.5*f2)
    dbTheta_dTht_m = - rSun**2*grid.sinTheta*br_m \
                     - 1j*m*bTor_m

    # I am overestimating the integral near the poles. Since the integrand
    #   is proportional to sinTheta, it goes to zero at the poles. Multiply 
    #   by 0.5 to get the area of the TRIANGLE insted of the RECTANGLE.
    #   Overestimation near the poles leads to instability in the diffusion
    #   term (calculated in rhs_diff).
    bTht_m[0] *= 0.5
    bTht_m[-1] *= 0.5

    return bTht_m, dbTheta_dTht_m

#def calcThetaComp(m, br_m, bTor_m):
#    
#    # This routine uses the trapezoidal rule to integrate
#    # IT DOES NOT WORK, CHECK
#
#    f1 = grid.sinTheta*br_m
#    f1[1:] = (f1[1:]+f1[:-1])*grid.dTheta/2.
#    f1[0] = grid.dTheta/4.*f1[0] # Boundary condition: f1 = 0 @ poles
#    f1[-1] *= 0.5 # Avoid numerical instability
#
#    f2 = bTor_m
#    f2[1:] = (f2[1:]+f2[:-1])*grid.dTheta/2.
#    f2[0] = grid.dTheta/4.*f2[0] # Boundary condition: f2 = 0 @ poles
#    f2[-1] *= 0.5 # Avoid numerical instability
#
#    bTht_m = - rSun**2*np.cumsum(f1) \
#             - 1j*m*np.cumsum(f2)
#
#    dbTheta_dTht_m = - rSun**2*grid.sinTheta*br_m \
#                     - 1j*m*bTor_m
#
#
#    return bTht_m, dbTheta_dTht_m


def rhs_advTheta(m, bTor_m):

    global mMax
    if True:#m<=mMax:
        rhs = -1./rSun*derivatives.deriv_theta(flows.uTheta*bTor_m)
    else:
        rhs = 0.

    return rhs


def rhs_shear(m, br_m, bTor_m, bTht_m, dbTht_dTheta_m):

    global mMax
    if True:#m<=mMax:
        
        term1 = grid.sinTheta*rSun**2*flows.omega_rSun*br_m
        term2 = flows.dOmega_dTheta_rNSSL*bTht_m
        term3 = flows.omega_rNSSL*dbTht_dTheta_m 

        rhs = term1 + term2 + term3
    else:
        rhs = 0.
              
    return rhs


def rhs_diff(m, br_m, bTor_m, bTht_m):
    
    global mMax
    if True:#m<=mMax:

        term1 = eta0/rSun**2*2.*grid.cosTheta/grid.sin3theta*1j*m*bTht_m
        term2 = -eta0/rSun**2/grid.sin2theta*bTor_m

        rhs = term1 + term2
    else:
        rhs = 0.

    return rhs


def advance(m, br_m, bTor_m, dtime):

    bTht_m, dbTheta_dTht_m = calcThetaComp(m, br_m, bTor_m)
    bTor_m_1 = bTor_m + rhs_shear(m, br_m, bTor_m, bTht_m, dbTheta_dTht_m)*dtime
    bTor_m_1 = bTor_m_1 + rhs_advTheta(m, bTor_m_1)*dtime

    bTor_m_1 = cranktor.crankStep(bTor_m_1, eta0, dtime, m) 
    bTht_m, dbTheta_dTht_m = calcThetaComp(m, br_m, bTor_m_1)
    bTor_m_1 = bTor_m_1 + rhs_diff(m, br_m, bTor_m_1, bTht_m)*dtime
    bTor_m_1 = bTor_m_1*np.exp(-1.*eta0/rSun**2/grid.sin2theta*m**2*dtime)

    if limitMagneticMemory:
        bTor_m_1 = bTor_m_1*np.exp(-dtime/decayTime)

    bTht_m, dbTheta_dTht_m = calcThetaComp(m, br_m, bTor_m_1)

    return bTor_m_1, bTht_m, dbTheta_dTht_m
