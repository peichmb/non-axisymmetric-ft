# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- timestep.py: timestep calculation 
# --- ----------------------------------------------------------------

import numpy as np
from params import *
import grid
import flows
tAcc = 0.

def getTimestep(n, nc):

    global tAcc

    if nSteps == 0:
        dtime = cfl*rSun*np.min(np.abs(grid.dTheta/flows.uTheta))
    else:
        dtime = timeBetMgs/float(nSteps)

    writeOut = False
    if nc%nStepsWrite == 0:
        writeOut = True

    endLoop = False
    #if tAcc + dtime >= timeBetMgs:
    if n%nSteps ==0:
        #dtime = timeBetMgs - tAcc
        #tAcc = 0.
        endLoop = True
    #else:
        #tAcc = tAcc + dtime

    return dtime, writeOut, endLoop
