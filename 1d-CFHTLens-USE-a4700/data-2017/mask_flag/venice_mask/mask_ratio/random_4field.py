import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

au = 1
def radec_range(field):
    if(field==1):
        dec1,dec2,ra1,ra2 = [(-11.244478)*au,(-3.6772671)*au , (30.17882)*au,(38.820652)*au]
    elif(field==2):
        dec1,dec2,ra1,ra2 = [-5.6952291*au,-0.952663*au ,132.0583*au,136.84*au]
    elif(field==3):
        dec1,dec2,ra1,ra2 = [51.19738*au,57.80508*au,208.55965*au,220.38112*au]
    else:
        dec1,dec2,ra1,ra2 = [-1.029174*au,4.6142831*au,329.97742*au,335.70685*au]
    return np.array([[ra1,ra2],[dec1,dec2]])

field = np.int32(np.sys.argv[1])
nrand = 3000000
bins = 500

# generate random positions in (ra,dec) space
xra = radec_range(field)
ra  = xra[0,0]+(xra[0,1]-xra[0,0])*np.random.ranf(nrand)
dec = xra[1,0]+(xra[1,1]-xra[1,0])*np.random.ranf(nrand)

# transfer (ra,dec) to (x,y)
arcs  = '1arcsec'
ff = fits.open('../W'+np.str(field)+'.16bit.'+arcs+'.reg2.fits')
w = wcs.WCS(ff[0].header,ff)

xy = w.wcs_world2pix(np.array([ra,dec]).T,1)
t=Table([  xy[:,0], xy[:,1], ra, dec ], names=('x','y','ra','dec'))
t.write('w'+np.str(field)+'_rand_xy.fits',format='fits')

ff.close()
