import numpy as np
import matplotlib.pylab as plt
import sys
field = np.int32(np.sys.argv[1])
zm    = np.sys.argv[2]
magl  = np.sys.argv[3]
rmax  = np.sys.argv[4]

#field,zm,magl,rmax = 1,'0.45','-21','0.025'#'0.0166667'
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
ratio = 0.1
nii,njj = 16,16 #for y and x
Grid = np.array([[8,9],[5,5],[7,7],[6,6]])
grid = Grid.copy()
grid[:,0],grid[:,1] = Grid[:,0]*nii,Grid[:,1]*njj
sg =np.array([(xr[0,1]-xr[0,0])/grid[field-1,0], (xr[1,1]-xr[1,0])/grid[field-1,1]])

#void
dirm = '/home/dfy/work/astro/1d-CFHTLens-USE-a4700/data-2017/mask_flag/maskmethod1_USE/w'+np.str(field)+'/'
ix,iy = 0,2
xv = np.genfromtxt(dirm+'w'+np.str(field)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos.dat')
xv[:,:2]=-xv[:,:2]

#grid void
print xv.shape
xvf = np.zeros((xv.shape[0],xv.shape[1]+1))
xvf[:,:-1] = xv
xvf[:,:2] = -xvf[:,:2]

#flag = [np.sum(nfield[:field-1])]
flag = [0]
nmean = xv.shape[0]*1./(Grid[field-1,0]*Grid[field-1,1])
#if(magl=='-21'):   nratio = 1
#if(magl=='-21.5'): nratio = 1.89
#if(magl=='-22'):   nratio = 2.77
#if(rmax=='0.0166667'): Nmean = 3000
#if(rmax=='0.025'):     Nmean = 500
#nmean = Nmean*nratio
nv = [0]
nij = [0]
def group_index(i,j,xr,xv,xvf,nv,Flag,flag,nij):
    for ii in np.arange(nii):
        for jj in np.arange(njj):
            #print nii*i+ii,njj*j+jj
            index1 = (xv[:,ix]>=xr[1,0]+sg[1]*(njj*j+jj))*(xv[:,ix]<xr[1,0]+sg[1]*(njj*j+jj+1))
            index2 = (xv[:,iy]>=xr[0,0]+sg[0]*(nii*i+ii))*(xv[:,iy]<xr[0,0]+sg[0]*(nii*i+ii+1))
            niijj = np.count_nonzero(index1*index2)
            #print niijj
            if(niijj>0) :
                nij[0] = nij[0] + niijj
                xvf[index1*index2,-1] = flag[0]
                if(nij[0]>=nmean*(1-ratio)):
                    #print nii*i+ii,njj*j+jj
                    #print nij[0]
                    Flag[nii*i+ii,njj*j+jj,0] = flag[0]
                    Flag[nii*i+ii,njj*j+jj,1] = nij[0]
                    flag[0] = flag[0]+1
                    nij[0] = 0
                nv[0] = nv[0]+niijj
                #print nv,niijj

Flag = np.zeros((grid[field-1,0],grid[field-1,1],2))
#for i in np.arange(1)+8:
for i in np.arange(Grid[field-1,0]):
    if(i%2==0):
    #for j in np.arange(1)+8: group_index(i,j,xr,xv,xvf,nv,flag,nij)
        for j in np.arange(Grid[field-1,1]): group_index(i,j,xr,xv,xvf,nv,Flag,flag,nij)
    if(i%2==1):
        for j in np.arange(Grid[field-1,1])[::-1]: group_index(i,j,xr,xv,xvf,nv,Flag,flag,nij)
print nv


np.savetxt(dirm+'w'+np.str(field)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos-flag-group3000-'+np.str(nii)+'-'+np.str(njj)+'.dat',xvf)

#plt
#plt.imshow(Flag[:,:,0])
#plt.imshow(Flag[:,:,1])

#plt.figure(figsize=(20,20))
#for i in np.arange(grid[field-1,0]):
#    for j in np.arange(grid[field-1,1]):
#        plt.axhline(xr[0,0]+sg[0]*j)
#        plt.axvline(xr[1,0]+sg[1]*i)
#        plt.text(xr[1,0]+sg[1]*i,xr[0,0]+sg[0]*j,Flag[i,j])
