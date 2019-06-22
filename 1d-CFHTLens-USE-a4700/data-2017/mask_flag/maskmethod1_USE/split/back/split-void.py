import numpy as np
import matplotlib.pylab as plt
import sys
nfield= np.array([72,25,49,25])
field = np.int32(np.sys.argv[1])
zm    = np.sys.argv[2]
magl  = np.sys.argv[3]
rmax  = np.sys.argv[4]

dirm = '/home/dfy/mw/dfy/CFHTLens/'
#grid
if(field == 1):
    xr = np.array([[-11.244478,-3.6772671],[30.17882,38.820652]])
if(field == 2):
    xr = np.array([[-5.6952291,-0.952663],[132.0583,136.84346]])
if(field == 3):
    xr = np.array([[51.19738,57.80508],[208.55965,220.38112]])
if(field == 4):
    xr = np.array([[-1.029174,4.6142831],[329.97742,335.70685]])

xr[:,0] = xr[:,0]-0.1
xr[:,1] = xr[:,1]+0.1
grid = np.array([[8,9],[5,5],[7,7],[6,6]])
sg =np.array([(xr[0,1]-xr[0,0])/grid[field-1,0], (xr[1,1]-xr[1,0])/grid[field-1,1]])

#void
dirm = '/home/dfy/mw/dfy/1d-CFHTLens-USE-a4700/mask_flag/maskmethod1_USE/w'+np.str(field)+'/'
ix,iy = 0,2
xv = np.genfromtxt(dirm+'w'+np.str(field)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos.dat')
xv[:,:2]=-xv[:,:2]

#grid void
print xv.shape
nv = 0.
xvf = np.zeros((xv.shape[0],xv.shape[1]+1))
xvf[:,:-1] = xv
xvf[:,:2] = -xvf[:,:2]
flag = np.sum(nfield[:field-1])
for i in np.arange(grid[field-1,0]):
    for j in np.arange(grid[field-1,1]):
        index1 = (xv[:,ix]>=xr[1,0]+sg[1]*j)*(xv[:,ix]<xr[1,0]+sg[1]*(j+1))
        index2 = (xv[:,iy]>=xr[0,0]+sg[0]*i)*(xv[:,iy]<xr[0,0]+sg[0]*(i+1))
        #print xr[1,0]+sg[1]*i, xr[1,0]+sg[1]*(i+1)
        nij = np.count_nonzero(index1*index2)
        if(nij>0) : 
            xvf[index1*index2,-1] = flag
            flag = flag+1
            #print nij
        nv = nv+nij
print nv
np.savetxt(dirm+'w'+np.str(field)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos-flag.dat',xvf)

#plt
#plt.figure(figsize=(8,8))
#step = 0.02#0.006139#0.2#0.006139
#nx0 = np.abs(np.ceil((xr[0,1]-xr[0,0])/step))
#ny0 = np.abs(np.ceil((xr[1,1]-xr[1,0])/step))
#h0,ybin0,xbin0=np.histogram2d(xv[:,iy],xv[:,ix],bins=[ny0,nx0],range=xr,normed=False)
#h0 = h0/(xbin0[1]-xbin0[0])/(ybin0[1]-ybin0[0])/(60**2.)
#plt.imshow(h0,interpolation='none',origin='low',extent=[xbin0[0], xbin0[-1], ybin0[0], ybin0[-1]])
#
#for i in np.arange(grid[field-1,0]):
#    plt.axhline(xr[0,0]+i*sg[0])
#for i in np.arange(grid[field-1,1]):
#    plt.axvline(xr[1,0]+i*sg[1])
