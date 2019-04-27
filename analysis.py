import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from sunpy.sun import carrington_rotation_number
from scipy.signal import argrelextrema


pi = np.pi
minSize = 100.
fname = 'totaldata_2.dat.cr.2'
pathSims = 'out/'

#pathSims = '/home/belda/toroidal/eqgrid/250.500/'
f = open(fname,'r')
lines = f.readlines()
f.close()
nPhi = 360
nTheta = 180
dTheta = 1.*pi/180.
m = np.fft.rfftfreq(nPhi)*nPhi

def deriv_theta(f):
    
    global nTheta
    f_theta = np.empty(nTheta, dtype='float')
    # Centered finite differences
    f_theta[1:nTheta-1] = (f[2:nTheta] - f[0:nTheta-2])/2./dTheta

    # Boundary conditions: f=0 @ poles
    f_theta[0] = (f[0] + f[1])/2./dTheta
    f_theta[nTheta-1] = -(f[nTheta-1] + f[nTheta-2])/2./dTheta

    return f_theta

firstc = 0
secondc = 0
thirdc = 0
fourthc = 0

for i in range(len(lines)):
    bits = lines[i].split()
    date = bits[0]
    year, month, day = date[:4], date[4:6], date[6:]
    date = year+'-'+month+'-'+day
    print date
    cr = int(carrington_rotation_number(date))
    if cr>1628:
        break

lines = lines[i:]

record = np.empty([len(lines),4])
i = 0
for line in lines:
    bits = line.split()
    date = bits[0]
    year, month, day = date[:4], date[4:6], date[6:]
    date = year+'-'+month+'-'+day
    cr = int(carrington_rotation_number(date))
    lat, lon, area = float(bits[2]), float(bits[3]), float(bits[4])
    record[i,:] = cr, lat, lon, area
    i += 1

del lines

bm0 = []
bmm = []
diff = []

for cr in range(1628,2007):

    if not cr in [1640,1641,1642,1643,1644,1645,1646,1647,1852,1853,1854,1855, 1856, 1857]:

        fname = 'm'+str(cr-2)+'f.bTor.fits'
        hl = fits.open(pathSims+fname)
        bPhi = hl[0].data
        hl.close()

        subRecord = record[np.where(record[:,0]==cr)]
        subRecord = subRecord[np.where(subRecord[:,3]>minSize)]

        lats = subRecord[:,1]
        lons = subRecord[:,2]
        areas = subRecord[:,3]
        
        w=np.where(np.abs(lats) <10)
        lats=lats[w]
        lons = lons[w]
        areas = areas[w]
        colats = 90.+lats

        for row in subRecord:
            coLat = int(row[1])+90
            lon = int(row[2])
            
            if np.abs(row[1]) < 10000.:
                    
                
                if lon == 360:
                    lon = 0
                bPhiStrip = bPhi[coLat, :]-np.mean(bPhi[coLat-5:coLat+5,:]+bPhi[179-coLat-5:179-coLat+5,:])/2.
                #bPhiStrip_fs = np.fft.rfft(bPhiStrip)
                # Get m = 0 component
                bPhiStrip_0 = np.mean(bPhiStrip)
                # Get the rest of m's 
                #bPhiStrip_m = np.zeros(bPhiStrip_fs.shape, dtype = 'complex')
                #bPhiStrip_m[1:] = bPhiStrip_fs[1:]
                #bPhiStrip_m = np.fft.irfft(bPhiStrip_m)
                bm0.append(bPhiStrip_0)
                bmm.append(bPhiStrip[lon]-bPhiStrip_0)
                print bPhiStrip_0, bPhiStrip[lon]
                if bm0[-1]>0 and bmm[-1]>0:
                    firstc = firstc+1
                elif bm0[-1]<0 and bmm[-1]>0:
                    secondc += 1
                elif bm0[-1]<0 and bmm[-1]<0:
                    thirdc += 1
                else:
                    fourthc += 1

            

bm0 = np.array(bm0)/1.e23
bmm = np.array(bmm)/1.e23

ap0 = bmm[np.where(bm0>0)]
an0 = bmm[np.where(bm0<0)]
print np.mean(ap0), np.mean(an0)

print firstc, secondc, thirdc, fourthc

eta = pathSims[-7:-5]
#eta = pathSims[-8:-5]
rflow = '5'
fig = plt.figure()
ax = fig.add_subplot(111)
title = r'$\eta = '+eta+'\,\mathrm{km^2 s^{-1}};\,u_r = '+rflow+'\,\mathrm{cm\,s^{-1}}.$'
plt.title(title, fontsize = 20)
plt.scatter(bm0, bmm, marker='x', color='green')
plt.axvline(0, color='black')
plt.axhline(0, color='black')
mm0 = np.max(np.abs(bm0))*1.2
mmm = np.max(np.abs(bmm))*1.2
plt.xlim([-mm0,mm0])
plt.ylim([-mmm,mmm])
plt.xlabel(r'$b_\phi^0\,\mathrm{[10^{23} Mx]}$', fontsize=17)
plt.ylabel(r'$b_\phi^m\,\mathrm{[10^{23} Mx]}$', fontsize=17)
ax.tick_params(labelsize=13)
plt.tight_layout()
#plt.savefig('comp.'+eta+'.'+rflow+'.png')
plt.show()

#bins = np.linspace(0,180,181)
#bincenters = np.linspace(0.5,179.5,180)
#hist, bins = np.histogram(diff,bins)
#plt.plot(bincenters,hist)
#plt.show()
#

#plt.show()
