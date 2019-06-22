import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

field = np.sys.argv[1]
zmean = np.sys.argv[2]
magl  = np.sys.argv[3]

ff = fits.open('w'+field+'random_xy.fits')
arcs  = '1arcsec'
ww = fits.open('../W'+np.str(field)+'.16bit.'+arcs+'.reg2.fits')
w = wcs.WCS(ww[0].header,ww)

index = ff[1].data['flag']>0
xy = w.wcs_pix2world(np.array([ff[1].data['x'][index],ff[1].data['y'][index]]).T,1)
xyz = np.zeros((xy.shape[0],3))
xyz[:,:2] = xy
xyz[:,2]  = np.ones(xy.shape[0])*np.float32(zmean)

np.savetxt('w'+field+'-z'+zmean+'-magl'+magl+'random.dat',xyz)
