#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fftw3.h"
#include "omp.h"

#include "func.h"
#include "struct.h"
#include "clusters.h"

int *thread_assign(int nthread,int n)
{
    int i,*randid,step,mod;

    randid = mycalloc(sizeof(int),nthread+1);
    step   = n/nthread;
    mod    = n%nthread;
    for(i=0;i<nthread;i++)
        randid[i] = step*i;
    //randid[nthread] = n-1;
    randid[nthread] = n;
    for(i=0;i<=nthread;i++)
        printf("%d ",randid[i]); printf("ngroup=%d\n",n);
    return randid;
}

void check_shear(SHEAR **shear,int n1,int n2)
{
    int i,j;
    for(i=0;i<n1;i++)
        for(j=0;j<n2;j++)
        {
            if(shear[i][j].gt!=0 || shear[i][j].gr!=0 || shear[i][j].wt!=0)
                printf("use of calloc is wrong!\n");
        }
}

SHEAR ****shear_thread_assign(int nzbin,int nthread,int ntr[2])
{
    int i,j;
    SHEAR ****shear;
    shear = mymalloc(sizeof(SHEAR ***)*nzbin);
    for(i=0;i<nzbin;i++)
    {
        shear[i] = mymalloc(sizeof(SHEAR **)*nthread);
        for(j=0;j<nthread;j++)
        {
            shear[i][j]=shear_calloc(ntr[0],ntr[1]);
            check_shear(shear[i][j],ntr[0],ntr[1]);
        }
    }
    return shear;
}

SHEARR ***shearr_thread_assign(int nzbin,int nthread,int ntr[2])
{
    int i;
    SHEARR ***shearr;
    shearr = mymalloc(sizeof(SHEARR **)*nzbin);
    for(i=0;i<nzbin;i++)
        shearr[i] = shearr_calloc(nthread,ntr[1]);
    return shearr;
}

void zero_shear_thread(SHEAR ****shear,SHEARR ***shearr,int nzbin,int ntr[2],int nthread)
{
    int i,j,k;
    for(k=0;k<nzbin;k++)
    {
        for(i=0;i<nthread;i++)
        {
            for(j=0;j<ntr[0];j++)
                memset(shear[k][i][j],0,sizeof(SHEAR)*ntr[1]);
            memset(shearr[k][i],0,sizeof(SHEARR)*ntr[1]);
        }
    }
}

void add_shear_thread(int nthread,SHEAR ***shear,int ntr[2])
{
    int i,j,k;
    for(i=1;i<nthread;i++)
        for(j=0;j<ntr[0];j++)
            for(k=0;k<ntr[1];k++)
            {
                shear[0][j][k].gt +=shear[i][j][k].gt;
                shear[0][j][k].gr +=shear[i][j][k].gr;
                shear[0][j][k].dr +=shear[i][j][k].dr;
                shear[0][j][k].angle +=shear[i][j][k].angle;
                shear[0][j][k].m +=shear[i][j][k].m;
                shear[0][j][k].err1 +=shear[i][j][k].err1;
                shear[0][j][k].err2 +=shear[i][j][k].err2;
                shear[0][j][k].nrr +=shear[i][j][k].nrr;
                shear[0][j][k].wt +=shear[i][j][k].wt;
            }

}

void zero_nfind_thread(int ***nfind,int n1,int ntr[2])
{
    int i,j,k;
    for(i=0;i<n1;i++)
        for(j=0;j<ntr[0];j++)
            for(k=0;k<ntr[1];k++)
            {
                nfind[i][j][k] = 0;
            }
}

void add_nfind_thread(int ****nfind,int nthread,int nzbinvoid,int ntr[2])
{
    int i,j,k,l;
    for(i=1;i<nthread;i++)
        for(l=0;l<nzbinvoid;l++)
            for(j=0;j<ntr[0];j++)
                for(k=0;k<ntr[1];k++)
                {
                    nfind[0][l][j][k] +=nfind[i][l][j][k];
                }
}

void add_nuse_thread(int nuse[],int nthread)
{
    int i;

    for(i=1;i<nthread;i++)
        nuse[0] +=nuse[i];
}
