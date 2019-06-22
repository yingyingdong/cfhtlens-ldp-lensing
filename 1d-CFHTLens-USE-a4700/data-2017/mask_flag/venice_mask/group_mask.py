#load astropy
import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits
from astropy.table import Table

field = '4'
#load (ra,dec) for galaxies in w1 shear catalogue
xx = np.genfromtxt('../w'+field+'_group_pos.dat')
radec = np.array((-xx[:,2],xx[:,1])).T
xy = np.zeros((xx.shape[0],2))

#load reg.fits and read its header
#arcs  = '5arcs'
arcs  = '1arcsec'
ff = fits.open('../W'+field+'.16bit.'+arcs+'.reg2.fits')
w = wcs.WCS(ff[0].header,ff)

#the 1-500000 galaxies in W1 catalogue
nn1,nn2 = 0,xx.shape[0]#500000
xy[nn1:nn2] = w.wcs_world2pix(radec[nn1:nn2],1)
#np.savetxt('w'+field+'_xy.cat',xy)
#t=Table([  xy[nn1:nn2,0], xy[nn1:nn2,1] ], names=('x','y'))
t=Table([  xy[nn1:nn2,0], xy[nn1:nn2,1], radec[nn1:nn2,0], radec[nn1:nn2,1] ], names=('x','y','ra','dec'))
t.write('w'+field+'_xy.fits',format='fits')

#use reg2.fits to judge wether galaxies are in the mask
#pix = np.int32(np.floor(xy[nn1:nn2]))
#mask  = ff[0].data
#mask_flag = mask[pix[:,0],pix[:,1]]
#print np.count_nonzero(mask_flag>1)*1./(nn2-nn1)

ff.close()
#np.savetxt('rand.cat',np.array(xy[:nn])[:,0,:])
