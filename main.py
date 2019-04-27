# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- David Martin-Belda, 2016
# --- ----------------------------------------------------------------

import numpy as np
from mpi4py import MPI
import aux
from params import *
import grid
import toroidal
import surface
import timestep
import output
import matplotlib.pyplot as plt

# Initialize MPI
comm = MPI.COMM_WORLD
nProc = comm.Get_size()
rank = comm.Get_rank()

# Set local matrix parameters
nFreqs = nPhi/2+1
chunkSizes, startPos, rowLims = aux.bufferTuples(nFreqs, nTheta, nProc)
m_local = grid.m[ rowLims[rank][0]:rowLims[rank][1] ]

# Time step
dtime = timeBetMgs/float(nSteps)

# Initial condition for br
#br = np.copy(br0)
br = np.zeros([nTheta, nPhi])
br_fs = aux.getMapFS(br)
br_fs_local = br_fs[ rowLims[rank][0]:rowLims[rank][1] , :]

# Initial condition for bTor
bTor = np.zeros([nTheta, nPhi])
bTor_fs = aux.getMapFS(bTor)
bTor_fs_local = bTor_fs[ rowLims[rank][0]:rowLims[rank][1] , :]

# Array to store bTht_fs
bTht_fs = np.zeros(bTor_fs.shape, dtype='complex')
bTht_fs_local = bTht_fs[ rowLims[rank][0]:rowLims[rank][1] , :]
dbTht_dTheta_fs = np.zeros(bTor_fs.shape, dtype='complex')
dbTht_dTheta_fs_local = dbTht_dTheta_fs[ rowLims[rank][0]:rowLims[rank][1] , :] 

# Arrays to store gathered maps
if rank == 0:
    br1_fs = np.empty(br_fs.shape, dtype='complex')
    if not solveOnlySFT:
        bTor1_fs = np.empty(bTor_fs.shape, dtype='complex')
        bTht1_fs = np.empty(bTht_fs.shape, dtype='complex')
        dbTht_dTheta1_fs = np.empty(bTht_fs.shape, dtype='complex')
else:
    br1_fs = None
    if not solveOnlySFT:
        bTor1_fs = None
        bTht1_fs = None
        dbTht_dTheta1_fs = None

if rank == 0:
    print
    print '-------- ---- ------- ----'
    print 'Toroidal flux mapping code'
    print '-------- ---- ------- ----'
    print
    print 'Running simulation on ' + str(nProc) + ' processes.'
    print
    print 'Surface meridional flow [cm/s]: ' + str(mflow_surf)
    print 'Return meridional flow [cm/s]: ' + str(mflow_return)
    print 'Diffusivity [cm**2/s]: ' + str(eta0)
    print 'Sources: '+sourcesPath
    print
    print '-------- ---- ------- ----'
    print

time = 0.
nWrite = 0



for elNumber in range(len(prefixList)):

    fname = prefixList[elNumber]
    currentNumber = int(fname)
    if rank == 0:
        print 'Evolving '+fname+' for '+str(timeToNextSources[elNumber])+' CRs...'

    if not solveOnlySFT:
        # Add br sources
        dbr = aux.readMap(sourcesPath+fname+'.dbr.fits')
        dbr_fs = aux.getMapFS(dbr)
        dbr_fs_local = dbr_fs[ rowLims[rank][0]:rowLims[rank][1] , :]
        br_fs_local += dbr_fs_local

        # Add bPhi sources
        dbPhi = aux.readMap(sourcesPath+fname+'.dbPhi.fits')
        dbPhi_fs = aux.getMapFS(dbPhi)
        dbPhi_fs_local = dbPhi_fs[ rowLims[rank][0]:rowLims[rank][1] , :]
        bTor_fs_local += dbPhi_fs_local
    else:
        br = aux.readMap(mgramsPath+fname+'.res.fits')
        br_fs = aux.getMapFS(br)
        br_fs_local = br_fs[ rowLims[rank][0]:rowLims[rank][1] , :]

    # Pack elements
    mPack = [ [m_local[i], br_fs_local[i], bTor_fs_local[i], bTht_fs_local[i], dbTht_dTheta_fs_local[i]] for i in range(len(m_local))]


    # Evolve each Fourier component according to the equations

    for nRots in range(timeToNextSources[elNumber]):

        endLoop = False
        n = 0
        while not endLoop:


            for item in mPack:

                writeOut = False
                nc = 0
                while not writeOut:#n<12000:#False:#not n<nSteps:#n<12000*12:#nSteps:#

                    m, br_m, bTor_m, bTht_m, dbTht_dTheta_m = item
                    nc += 1
                    # Advance toroidal field
                    if not solveOnlySFT:#False:#
                        bTor_m[:], bTht_m[:], dbTht_dTheta_m[:] = toroidal.advance(m, br_m, bTor_m, dtime)
                    # Advance surface field
                    br_m[:] = surface.advance(m, br_m, dtime)
                    if nc%nStepsWrite == 0:
                        writeOut = True



            n += nStepsWrite

            if n%nSteps == 0:
                endLoop = True

            #if rank == 0:
            #    print m, n, nc, nStepsWrite, writeOut, endLoop

            comm.Gatherv(br_fs_local, [br1_fs, chunkSizes, startPos, MPI.DOUBLE_COMPLEX])
            if not solveOnlySFT:
                comm.Gatherv(bTor_fs_local, [bTor1_fs, chunkSizes, startPos, MPI.DOUBLE_COMPLEX])
                #comm.Gatherv(bTht_fs_local, [bTht1_fs, chunkSizes, startPos, MPI.DOUBLE_COMPLEX])
            # Write output
            if rank == 0:
                if not solveOnlySFT:
                    if TEST:
                        fileName = outPath+fname+'.'+str(nWrite).rjust(5,'0')
                    else:
                        fileName = outPath+str(currentNumber+nRots)+'.'+str(n).rjust(len(str(nSteps)),'0')
                        print n, nSteps, fileName
                    br1 = aux.getMapRS(br1_fs)
                    bTor1 = aux.getMapRS(bTor1_fs)
                    #bTht1 = aux.getMapRS(bTht1_fs)
                    output.saveMap(br1, fileName+'.br.fits')
                    output.saveMap(bTor1, fileName+'.bTor.fits')
                    #output.saveMap(bTht1, fileName+'.bTht.fits')
                    print '-- Output written at time [years]: '+str(round(time/24./3600./365.,3))
                else:
                    if nRots == timeToNextSources[elNumber]-1:
                        br1 = aux.getMapRS(br1_fs)
                        output.saveMap(br1, outPath+fname+'.res.evolved.fits')
                        print 'Wrote '+fname+'.res.evolved.fits'

            nWrite += 1



        time = time + timeBetMgs

        
if rank == 0:
    print '---'
    print 'Simulation ended at time = ' + str(round(time/365./24./3600.,3)) + ' years.'
    print 'Output stored in '+outPath
    print
