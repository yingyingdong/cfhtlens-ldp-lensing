#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fftw3.h"
#include "func.h"
#include "struct.h"
#include "Romberg.h"
#include "power.h"
#include "clusters.h"

double *zbin_lg(int nzbinlg,double zmin,double zmax,double *zstep)
{
    int i;
    double *zbinlg,step;

    zbinlg = mymalloc(sizeof(double)*nzbinlg);
    step = (log10(1+zmax)-log10(1+zmin))/nzbinlg;
    for(i=0;i<nzbinlg;i++)
        zbinlg[i] = pow(10.,(i+0.5)*step+log10(1+zmin))-1;
    *zstep = step;
    return zbinlg;
}

double gain_zbin_lg(double *DAzlg,double z,double zmin,double zstep)
{
    int nz;
    nz = floor( (log10(1+z)-log10(1+zmin))/zstep );
    return DAzlg[nz];
}

double DA_z_ww(double z)
{
    return DH*qromb(angu_dia_d_ww,0,z)/(1.+z);
}

double D_COM_ww(double z)
{
    return DH*qromb(angu_dia_d_ww,0,z);
}


double *DA_zlg_ww(int nzbinlg,double *zbinlg)
{
    int i,j;
    double *DAz;

    DAz = mymalloc(sizeof(double)*nzbinlg);
    for(i=0;i<nzbinlg;i++)
    {
        DAz[i] = DA_z_ww(zbinlg[i]);
        //printf("%f %f\n",zbinlg[i],DAz[i]);
    }
    return DAz;
}


//double DA_z12(double result1,double result2,double z1)
//{
//    return 1./(1+z2)*(result2-result1);
//}

double sigma(double z1,double z2,double zmin,double zstep,double DA1,double *DAz)
{
    double DA2,DA12;
    int nz1,nz2;

    //nz1 = floor( (log10(1+z1)-log10(1+zmin))/zstep );//floor((z2-zmin)/zstep);
    nz2 = floor( (log10(1+z2)-log10(1+zmin))/zstep );//floor((z2-zmin)/zstep);
    //DA1 = DAz[nz1];
    DA2 = DAz[nz2];
    DA12 = 1./(1+z2)*(DA2*(1+z2)-DA1*(1+z1));
    //printf("DA1=%f,DA2=%f,DA12=%f,sigma=%f\n",DA1,DA2,DA12,cM*cM/(4*pi*GM)*(DA1/1000));
    return cM*cM/(4*pi*GM)*(DA2/1000.)/((DA1/1000.)*(DA12/1000.));  //Mpc
}
