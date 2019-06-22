import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits

nxy = [[1266,1441],[806,813],[1109,1959],[952,966]]

nfield = np.int32(np.sys.argv[1])
field = np.str(nfield)
xx = fits.open('w'+field+'_xy_mask.fits')
radec = np.array((-xx[1].data['ra'],xx[1].data['dec'],xx[1].data['flag'])).T
flag = radec[:,2].reshape(nxy[nfield-1][0],nxy[nfield-1][1])
#plt.imshow(flag)
np.savetxt('w'+field+'_xy_mask.dat',flag,fmt='%1d')
