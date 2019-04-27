# --- ----------------------------------------------------------------
# --- Code to calculate the latitudinal and longitudinal distributions
# ---   of toroidal flux from magnetograms.
# --- output.py: contains storing and plotting routines 
# --- ----------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from params import *
import grid
import cPickle as cP


def mkfigs(br, bTor, bTht, dbTht):

    plt.figure()
    plt.title('Radial surface field')
    mm = np.max(np.abs(br))
    plt.imshow(br, cmap = 'RdBu_r', vmin=-mm, vmax=mm)
    plt.colorbar(orientation = 'horizontal')

    plt.figure()
    plt.title('Integrated toroidal field')
    mm = np.max(np.abs(bTor))
    plt.imshow(bTor, cmap = 'RdBu_r',vmin=-mm, vmax=mm)
    plt.colorbar(orientation = 'horizontal')

    plt.figure()
    plt.title('bTheta')
    mm = np.max(np.abs(bTht))
    plt.imshow(bTht, cmap = 'RdBu_r',vmin=-mm, vmax=mm)
    plt.colorbar(orientation = 'horizontal')

    plt.figure()
    plt.title('dbTheta_dTheta')
    mm = np.max(np.abs(dbTht))
    plt.imshow(dbTht, cmap = 'RdBu_r',vmin=-mm, vmax=mm)
    plt.colorbar(orientation = 'horizontal')

    plt.show()


def saveMap(b, fname):

    fits.writeto(fname, b)

