import numpy as np
import matplotlib.pylab as plt
import sys
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
nx,ny = 100,100
ratio = 0.1
nii,njj = 2,2 #for y and x
grid = np.array([[nx,ny]])
Grid = np.array([[nx/nii,ny/njj]])
sg =np.array([(xr[0,1]-xr[0,0])/grid[field-1,0], (xr[1,1]-xr[1,0])/grid[field-1,1]])
nfield= np.array([grid[field-1,0]*grid[field-1,1]])
print nfield

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
nmean = xv.shape[0]*1./(Grid[field-1,0]*Grid[field-1,1])
nij = 0
def group_index(i,j,xr,xv,xvf,nv,flag):
    for ii in np.arange(nii)
        for jj in np.arange(njj):
            index1 = (xv[:,ix]>=xr[1,0]+sg[1]*(njj*j+jj))*(xv[:,ix]<xr[1,0]+sg[1]*(njj*j+jj+1))
            index2 = (xv[:,iy]>=xr[0,0]+sg[0]*(nii*i+ii))*(xv[:,iy]<xr[0,0]+sg[0]*(nii*i+ii+1))
            niijj = np.count_nonzero(index1*index2)
            if(niijj>0) :
                nij = nij+niijj
                xvf[index1*index2,-1] = flag
                if(nij>=nmean*(1-ratio)):
                    nij = 0
                    flag = flag+1
                    Flag[nii*i+ii,njj*j+jj] = flag
                nv = nv+niijj

Flag = np.zeros((grid[field-1,0],grid[field-1,1]))
for i in np.arange(Grid[field-1,0]):
    if(j/2==0):
        for j in np.arange(Grid[field-1,1]): group_index(i,j,xr,xv,xvf,nv,flag)
    if(j/2==1):
        for j in np.arange(Grid[field-1,1])[::-1]: group_index(i,j,xr,xv,xvf,nv,flag)
print nv


#np.savetxt(dirm+'w'+np.str(field)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos-flag'+np.str(nfield[field-1])+'.dat',xvf)

#plt
plt.figure(figsize=(8,8))
for i in np.arange(grid[field-1,0]):
    for j in np.arange(grid[field-1,1]):
        plt.axhline(grid[field-1,0])
        plt.axvline(grid[field-1,1])
        plt.txt(grid[field-1,0]*(1+0.5),grid[field-1,1]*(1+0.5),flag[i,j])
