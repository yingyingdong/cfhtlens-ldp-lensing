x = np.genfromtxt('w1_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean0.45-rmax0.0166667_magl-21-maskratio0.1-pos-flag-group3000-16-16.dat')
ima,edgex,edgey=np.histogram2d(x[:,0],x[:,2],bins=(20,20),normed=True)
nx,ny = edgex.shape[0]-1,edgey.shape[0]-1

flag = np.zeros((ny,nx))
for i in np.arange(nx):
    for j in np.arange(ny):
        index = (x[:,0]>edgex[i])*(x[:,0]<edgex[i+1])*(x[:,2]>edgey[j])*(x[:,2]<edgey[j+1])
        flag[i,j] = np.mean(x[index,-1])

plt.imshow(flag,interpolation='none',origin='low',extent=(edgey[0],edgey[-1],edgex[0],edgex[-1]))
