import numpy as np
import astropy
from astropy import wcs
from astropy.io import fits
from astropy.table import Table
import matplotlib.pylab as plt

field = np.sys.argv[1]
x = fits.open('w'+field+'_rand_xy_mask.fits')
index = x[1].data['flag']==1
bins = 500
r1,r2 = np.count_nonzero(index),x[1].data['flag'].shape[0]
print(r1*1./r2)
plt.hist2d(x[1].data['ra'][index],x[1].data['dec'][index],bins=(bins,bins))
