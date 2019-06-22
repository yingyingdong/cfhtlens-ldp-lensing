import numpy as np
import matplotlib.pylab as plt
import sys
zm    = np.sys.argv[1]
magl  = np.sys.argv[2]
rmax  = np.sys.argv[3]

Nmean = 3000
nii,njj = 16,16
idf  = -1
flag = 0

dirm = '/home/dfy/work/astro/1d-CFHTLens-USE-a4700/data-2017/mask_flag/maskmethod1_USE/'
for ifield in np.arange(4)+1:
    xvf = np.genfromtxt(dirm+'w'+np.str(ifield)+'/'+'w'+np.str(ifield)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos-flag-group'+np.str(Nmean)+'-'+np.str(nii)+'-'+np.str(njj)+'.dat')
    nflag = np.max(xvf[:,idf])
    xvf[:,idf] = xvf[:,idf] + flag
    flag  = flag + nflag
    print nflag, flag
    np.savetxt(dirm+'w'+np.str(ifield)+'/'+'w'+np.str(ifield)+'_maskmethod1_nzbin1-step0.006139-zstep0.2-zmean'+zm+'-rmax'+rmax+'_magl'+magl+'-maskratio0.1-pos-flag-group'+np.str(Nmean)+'-'+np.str(nii)+'-'+np.str(njj)+'final.dat',xvf)
