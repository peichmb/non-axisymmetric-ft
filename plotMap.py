import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from sys import argv

fname = argv[1]

hl = fits.open(fname)
b = hl[0].data
hl.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
mm = np.max(np.abs(b))
plt.imshow(b, vmin=-mm, vmax=mm, cmap='RdBu_r')
plt.xlabel('Longitude [deg]')
plt.ylabel('Colatitude [deg]')
plt.colorbar(orientation='horizontal')
plt.show()

