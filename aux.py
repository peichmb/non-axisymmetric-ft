# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- aux.py: This module contains auxiliary functions, for tasks such
# ---   as reading the magnetograms or calculating the parameters for
# ---   collective operations.
# --- ----------------------------------------------------------------

import numpy as np
import astropy.io.fits as fits
import cPickle as cP


def readMap(fname):
    
    hl = fits.open(fname)
    map = hl[0].data
    hl.close()

    return map


def getMapFS(map):
    nTheta, nPhi = map.shape
    nFreqs = nPhi/2+1
    map_fs = np.empty([nFreqs,nTheta], dtype='complex')
    for j in range(nTheta):
        map_fs[:,j] = np.fft.rfft(map[j,:])

    return map_fs
    

def getMapRS(map_fs):

    nFreqs, nTheta = map_fs.shape
    nPhi = (nFreqs-1)*2
    map = np.empty([nTheta,nPhi])
    for i in range(nTheta):
        map[i,:] = np.fft.irfft(map_fs[:,i])

    return map


def bufferTuples(nRows,nCols,nProc):

    hChunks0 = nRows/nProc
    nExtra = nRows%nProc
    
    sizes = [hChunks0 for i in range(nProc)]
    sizes[:nExtra] = [hChunks0+1 for i in range(nExtra)]

    ll = [0]
    uu = []
    for i in sizes[:-1]:
        ll.append(ll[-1]+i)
    uu = ll[1:]
    uu.append(nRows)
    rowLims = []
    for i in range(nProc):
        rowLims.append((ll[i],uu[i]))
    rowLims = tuple(rowLims)

    sizes = [i*nCols for i in sizes]
    sizes = tuple(sizes)

    startPos = [0]
    for i in sizes[:-1]:
        startPos.append(startPos[-1]+i)
    startPos = tuple(startPos)
    
    return sizes, startPos, rowLims


def loadSrc(fname):

    f = fits.open(fname)
    b = np.float64(f[0].data)
    f.close()
    
    return b
