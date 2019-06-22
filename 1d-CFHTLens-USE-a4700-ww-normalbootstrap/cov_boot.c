#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fftw3.h"
#include "omp.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include "func.h"
#include "struct.h"
#include "Romberg.h"
#include "power.h"
#include "clusters.h"
#include "get_z_dal.h"
#include "thread_assign.h"
#include "output.h"

#define nr 20

void return_fname(char *fname,double magl,double zmean,double rcut,double dz,int mul,int iboot);
void load_boot_err(char *fname,double **xx);
void cal_cov(double ***xx,double **cov);
void out_cov(double **cov,double magl,double rcut,double dz,double zmean,int mul);

int main(int argc,char *argv)
{

    int i,j;
    double **cov,***xx; 
    double magl=-21.5,zmean=0.45,rcut=1,dz=0.2,mul=160000;
    char fname[1024];

    cov = D_malloc2(nr,nr);
    xx  = mymalloc(sizeof(double **)*nboot);
    for(i=0;i<nboot;i++)
        xx[i] = D_malloc2(nr,10);

    for(i=0;i<nboot;i++)
    //for(i=0;i<1;i++)
    {
        return_fname(fname,magl,zmean,rcut,dz,mul,i);
        load_boot_err(fname,xx[i]);
    }
    cal_cov(xx,cov);
    out_cov(cov,magl,rcut,dz,zmean,mul); 
    return 0;
}

void out_cov(double **cov,double magl,double rcut,double dz,double zmean,int mul)
{
    int i,j;
    char ofile[1024];
    FILE *fp;
    sprintf(ofile,"data-2017/figure/cov/4jzfield-ww%lg-magl%lg-rcut%lg-dz%lg-zmean%lg-mul%d",
            ww,magl,rcut,dz,zmean,mul);
    myfopen(fp,ofile,"w");
    for(i=0;i<nr;i++)
    {
        for(j=0;j<nr;j++)
            fprintf(fp,"%lg ",cov[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}


void return_fname(char *fname,double magl,double zmean,double rcut,double dz,int mul,int iboot)
{
    sprintf(fname,"data-2017/onetime-sigma-pz/4jzfield-ww%lg-void_mag%lg_zmean%lg-rcut%lg-nbin1-zstep%lg-step0.006139-zmean%lg-1-20-shear-cali2-lgadr-COM-mul%dcra-0.58-2.5-wtsigma-pz0.2omp-iboot%d.dat",ww,magl,zmean,rcut,dz,zmean,mul,iboot);
    printf("%s\n",fname);
}

void load_boot_err(char *fname,double **xx)
{
    int i,j;
    FILE *fp;
    myfopen(fp,fname,"r");
    for(i=0;i<nr;i++)
        for(j=0;j<10;j++)
            fscanf(fp,"%lg",&xx[i][j]);
    for(j=0;j<10;j++) printf("%lg ",xx[nr-1][j]); printf("\n");
    fclose(fp);
}

void cal_cov(double ***xx,double **cov)
{
    int i,j,k,ie=0;
    double *xmean;
    xmean = mymalloc(sizeof(double)*nr);

    for(i=0;i<nr;i++)
    {
        xmean[i] = 0;
        for(j=0;j<nboot;j++)
            xmean[i] +=xx[j][i][ie]/nboot;
    }

    for(i=0;i<nr;i++)
        for(j=0;j<nr;j++)
        {
            for(k=0;k<nboot;k++)
                cov[i][j] +=(nboot-1.)/nboot*(xx[k][i][ie]-xmean[i])*(xx[k][j][ie]-xmean[j]);
        }
}
