import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits

field = np.sys.argv[1]
arcs  = '5arcs'
#arcs  = '1arcsec'
fname = ('../W'+field+'.16bit.'+arcs+'.reg2.fits')
#fname = ('../W'+field+'.16bit.'+arcs+'.reg2.fits')
data, header = fits.getdata(fname, header=True)
index1 = data>1
data_w = data.copy()
data_w[index1] = 0
index2 = data<=1 
data_w[index2] = 1

#overwrite:
fits.writeto('output_w'+field+'.'+arcs+'.fits', data_w, header, clobber=True)
