import numpy as np
from sunpy import sun
import matplotlib.pyplot as plt
from astropy.io import fits
from sys import argv

fname = 'totaldata_2.dat.cr'
pathSims = 'out/'
pathMgs = 'resampled/'
f = open(fname,'r')
lines = f.readlines()
f.close()

pi = np.pi
nTheta = 180
nPhi = 360
dTheta = 1.*pi/180.

def deriv_theta(f):
    
    global nTheta
    f_theta = np.empty(nTheta, dtype='float')
    # Centered finite differences
    f_theta[1:nTheta-1] = (f[2:nTheta] - f[0:nTheta-2])/2./dTheta

    # Boundary conditions: f=0 @ poles
    f_theta[0] = (f[0] + f[1])/2./dTheta
    f_theta[nTheta-1] = -(f[nTheta-1] + f[nTheta-2])/2./dTheta

    return f_theta

for i in range(len(lines)):
    bits = lines[i].split()
    date = bits[0]
    year, month, day = date[:4], date[4:6], date[6:]
    date = year+'-'+month+'-'+day
    print date
    cr = int(sun.carrington_rotation_number(date))
    if cr>1874:
        break

lines = lines[i:]

crCurrent = 1875
colats = []
lons = []
ss = []

for i in range(len(lines)):

    bits = lines[i].split()
    date = bits[0]
    year, month, day = date[:4], date[4:6], date[6:]
    date = year+'-'+month+'-'+day
    print date
    cr = int(sun.carrington_rotation_number(date))
    colat = 90.+float(bits[2])
    lon = float(bits[3])
    s = float(bits[4])

    if cr in range(1640,1647) or cr in range(1853,1858):
        print 'Ignoring CR ',cr
        crCurrent = cr+1
    elif cr == crCurrent:
        if s>50:#
            colats.append(colat)
            lons.append(lon)
            ss.append(s)
            print cr, lon, colat
    else:
        cr = crCurrent
        print '-----' 
        crNext = str(cr+1)
        crPrev = str(cr-1)
        crPrevPhi = str(cr-2)
        cr = str(cr)
        print crCurrent, crNext, cr, crPrev

        hl = fits.open(pathMgs+'m'+crPrev+'f.res.fits')
        brPrev = hl[0].data
        hl.close()
        hl = fits.open(pathMgs+'m'+cr+'f.res.fits')
        br = hl[0].data
        hl.close()
        hl = fits.open(pathMgs+'m'+crNext+'f.res.fits')
        brNext = hl[0].data
        hl.close()
        hl = fits.open(pathSims+'m'+crPrevPhi+'f.bTor.fits')
        bPhiPrev = hl[0].data
        hl.close()

        # Remove axisymmetric component
        i1 = np.empty(br.shape)
        for nCol in range(nPhi):
            i1[:,nCol] = np.mean(bPhiPrev,1)
        i1r = np.empty(br.shape)
        for nRow in range(nTheta):
            i1r[nRow,:] = i1[nTheta-1-nRow,:]
        bPhiPrev = bPhiPrev - (i1+i1r)

        #plt.imshow(i1)
        #plt.figure()
        #plt.imshow(i1r)
        #plt.figure()
        #plt.imshow(i1+i1r)
        #plt.show()
        #exit()

        if True:
            bPhiPrev_der = np.empty(bPhiPrev.shape)
            for j in range(nPhi):
                bPhiPrev_der[:,j] = deriv_theta(bPhiPrev[:,j])
            #bPhiPrev = bPhiPrev_der


        fig = plt.figure(figsize=(5,9))
        ax = fig.add_subplot(311)
        plt.title('CR '+crPrev)
        #plt.xlabel('Longitude[deg]')
        plt.ylabel('Colatitude[deg]')
        mm = np.max(np.abs(brPrev))/4.
        plt.imshow(brPrev, cmap='RdBu_r', vmin = -mm, vmax = mm, origin='lower')
        plt.colorbar()
        plt.scatter([lons],[colats], color='green', marker='x')#, s = ss)
        plt.xlim([0,360])
        plt.ylim([0,180])
        #fig = plt.figure()
        ax = fig.add_subplot(312)
        plt.title('bPhi; CR '+crPrev)
        #plt.xlabel('Longitude[deg]')
        plt.ylabel('Colatitude[deg]')
        mm = 5.e23 
        plt.imshow(bPhiPrev, cmap='RdBu_r', origin='lower', vmin = -mm, vmax = mm)#, origin='lower')
        plt.colorbar()
        plt.scatter([lons],[colats], color='green', marker='x')
        plt.xlim([0,360])
        plt.ylim([0,180])
        #ax = fig.add_subplot(223)
        #plt.title('der bPhi; CR '+crPrev)
        #plt.xlabel('Longitude[deg]')
        #plt.ylabel('Colatitude[deg]')
        #mm = 3.e23 
        #mm = np.max(np.abs(bPhiPrev_der))
        #plt.imshow(bPhiPrev_der, cmap='RdBu_r', origin='lower', vmin = -mm, vmax = mm)#, origin='lower')
        #plt.colorbar()
        #plt.scatter([lons],[colats], color='green')
        #plt.xlim([0,360])
        #plt.ylim([0,180])
        ax = fig.add_subplot(313)
        plt.title('CR '+cr)
        plt.xlabel('Longitude[deg]')
        plt.ylabel('Colatitude[deg]')
        mm = np.max(np.abs(br))/4.
        plt.imshow(br, cmap='RdBu_r', vmin = -mm, vmax = mm, origin='lower')
        plt.colorbar()
        plt.scatter([lons],[colats], color='green', marker='x')
        plt.xlim([0,360])
        plt.ylim([0,180])
        #plt.tight_layout()
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #plt.title('CR '+crNext)
        #plt.xlabel('Longitude[deg]')
        #mm = np.max(np.abs(brNext))/4.
        #plt.ylabel('Colatitude[deg]')
        #plt.imshow(brNext, cmap='RdBu_r', vmin = -mm, vmax = mm, origin='lower')
        #plt.scatter([lons],[colats], color='green')
        #plt.tight_layout()
        #fig = plt.figure()
        plt.tight_layout()
        plt.savefig('frames/'+str(crCurrent)+'.png')
        #plt.show()

        #plt.figure()
        #plt.title('bPhi; CR '+crPrev)
        #plt.xlabel('Longitude[deg]')
        #plt.ylabel('Colatitude[deg]')
        #mm = 5.e23 
        #plt.imshow(bPhiPrev, cmap='RdBu_r', origin='lower', vmin = -mm, vmax = mm)#, origin='lower')
        #plt.colorbar(orientation = 'horizontal')
        #plt.scatter([lons],[colats], color='green')
        #plt.xlim([0,360])
        #plt.ylim([35,145])
        #plt.show()
        #raw_input()
        
        crCurrent += 1
        colats = []
        lons = []
        ss = []

