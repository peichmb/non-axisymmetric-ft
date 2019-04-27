# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- params.py: contains the parameters of the problem.
# --- ----------------------------------------------------------------

import aux
import numpy as np


# Constants and conversion factors
# --------------------------------

TEST = False

pi = np.pi
rSun = 695.8e8 # [cm]

days2secs = 24.*3600.
secs2days = 1./days2secs

years2secs = 365.*days2secs
secs2years = 1./years2secs

# What problem do I want to solve?
solveOnlySFT = False


# Problem parameters
# ------------------

eta0 = 500.e10 # Diffusion coefficient [cm**2/s]
mflow_return = 600. # Magnitude of return meridional flow [cm/s]

limitMagneticMemory = False
decayTime = 2.5 # Magnetic memory of the system
decayTime *= years2secs

eta_surf = 250.e10 # Surface turbulent diffusion
mflow_surf = 1100. # Magnitude of surface meridional flow

dRotModel = 2 # Differential rotation model (see flows.py)
timeBetMgs = 27.2753*days2secs # Time between magnetograms

nSteps = 1000
nStepsWrite = 40

# -- Get nPhi and nTheta from a sample magnetogram
map = aux.readMap('resampled/1625.res.fits')
nTheta, nPhi = map.shape
nx = nTheta
del map

# -- List of input magnetograms
nStart = 1625 # First magnetogram
nEnd = 2126 # Last magnetogram
prefixList = [str(i) for i in range(nStart,nEnd+1)]
# Remove missing magnetograms
for i in range(1640,1645):
    prefixList.pop(prefixList.index(str(i)))
prefixList.pop(prefixList.index('1854'))
prefixList.pop(prefixList.index('2015'))
prefixList.pop(prefixList.index('2016'))
prefixList.pop(prefixList.index('2040'))
prefixList.pop(prefixList.index('2041'))

# Time to next injection of sources in carrington rotations
timeToNextSources = [1 for i in prefixList]
timeToNextSources[prefixList.index('1639')] = 6
timeToNextSources[prefixList.index('1853')] = 2
timeToNextSources[prefixList.index('2014')] = 3
timeToNextSources[prefixList.index('2039')] = 3

# -- Paths
etaStr = str(int(eta0/1.e10)).rjust(4,'0')
mflowStr = str(int(mflow_return)).rjust(4,'0')
if limitMagneticMemory:
    decayStr = str(int(round(decayTime*secs2years*10))).rjust(4,'0')
else:
    decayStr = 'xxxx'
velModelStr = str(dRotModel)

outPath = '/scratch/belda/toroidal/sources_hathaway/new_mgs_fine2/'
outPath += velModelStr+'.'+etaStr+'.'+mflowStr+'.'+decayStr+'/'

# ---- full run
sourcesPath = 'sources_hathaway/'
# ---- only sft
mgramsPath = 'resampled/'

# TEST
# ----
if TEST:
    prefixList = ['test1']
    timeToNextSources = [13*20] # Integer number in units of timeBetMgs
    sourcesPath = 'sources_test/'
    outPath = '/scratch/belda/toroidal/tests/out_test1/'

