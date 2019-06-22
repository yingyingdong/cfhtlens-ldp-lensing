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

SHEAR ****shear_thread_assign_nboot(int nthread,int ntr[2])
{
    int j,k;
    SHEAR ****shear;
    shear = mymalloc(sizeof(SHEAR ***)*nthread);
    
        for(j=0;j<nthread;j++)
        {
            shear[j] = mymalloc(sizeof(SHEAR **)*(nboot+1));
            for(k=0;k<=nboot;k++)
            {
                shear[j][k] = shear_calloc(ntr[0],ntr[1]);
                check_shear(shear[j][k],ntr[0],ntr[1]);
            }
        }
    return shear;
}

SHEARR ***shearr_thread_assign_nboot(int nthread,int ntr[2])
{
    int i,j;
    SHEARR ***shearr;

        shearr = mymalloc(sizeof(SHEARR **)*nthread);
        for(i=0;i<nthread;i++)
            shearr[i] = shearr_calloc((nboot+1),ntr[1]);
    return shearr;
}

void zero_shear_thread_nboot(SHEAR ****shear,SHEARR ***shearr,int **nuse,int ntr[2],int nthread)
{
    int i,j,iboot;
        for(i=0;i<nthread;i++)
            for(iboot=0;iboot<=nboot;iboot++)
            {
                for(j=0;j<ntr[0];j++)
                    memset(shear[i][iboot][j],0,sizeof(SHEAR)*ntr[1]);
                memset(shearr[i][iboot],0,sizeof(SHEARR)*ntr[1]);
                nuse[i][iboot] = 0 ;
                //printf("shear[0][0][0][0].wt=%lg\n",shear[0][0][0][0].wt); fflush(stdout);
            }
}

void zero_shear_ntr(SHEAR **shear,int ntr[2])
{
    int i,j;
    for(i=0;i<ntr[0];i++)
        memset(shear[i],0,sizeof(SHEAR)*ntr[1]);
}

void add_shear_thread(int nthread,SHEAR ****shear,int ntr[2],int iboot)
{
    int i,j,k;
    for(i=1;i<nthread;i++)
        for(j=0;j<ntr[0];j++)
            for(k=0;k<ntr[1];k++)
            {
                shear[0][iboot][j][k].gt    +=shear[i][iboot][j][k].gt;
                shear[0][iboot][j][k].gr    +=shear[i][iboot][j][k].gr;
                shear[0][iboot][j][k].dr    +=shear[i][iboot][j][k].dr;
                shear[0][iboot][j][k].angle +=shear[i][iboot][j][k].angle;
                shear[0][iboot][j][k].m     +=shear[i][iboot][j][k].m;
                shear[0][iboot][j][k].err1  +=shear[i][iboot][j][k].err1;
                shear[0][iboot][j][k].err2  +=shear[i][iboot][j][k].err2;
                shear[0][iboot][j][k].nrr   +=shear[i][iboot][j][k].nrr;
                shear[0][iboot][j][k].wt    +=shear[i][iboot][j][k].wt;
            }

}

void add_shear_to_iboot(SHEAR **tmp_shear,SHEAR ***shear,FLAG_GROUP *flag_group,int *nuse,int ntr[2],int index)
{
    int iboot,j,k,nn;

    for(iboot=0;iboot<=nboot;iboot++)
        if(flag_group[iboot].flag[index]>0)
        {
            if(iboot<nboot)       nn = flag_group[iboot].flag[index];
            else if(iboot==nboot) nn = 1;
            for(j=0;j<ntr[0];j++)
                for(k=0;k<ntr[1];k++)
                {
                    shear[iboot][j][k].gt    +=nn*tmp_shear[j][k].gt;
                    shear[iboot][j][k].gr    +=nn*tmp_shear[j][k].gr;
                    shear[iboot][j][k].dr    +=nn*tmp_shear[j][k].dr;
                    shear[iboot][j][k].angle +=nn*tmp_shear[j][k].angle;
                    shear[iboot][j][k].m     +=nn*tmp_shear[j][k].m;
                    shear[iboot][j][k].err1  +=nn*tmp_shear[j][k].err1;
                    shear[iboot][j][k].err2  +=nn*tmp_shear[j][k].err2;
                    shear[iboot][j][k].nrr   +=nn*tmp_shear[j][k].nrr;
                    shear[iboot][j][k].wt    +=nn*tmp_shear[j][k].wt;
                }
            nuse[iboot] ++;
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

void add_nuse_thread(int **nuse,int nthread,int iboot)
{
    int i;

    for(i=1;i<nthread;i++)
        nuse[0][iboot] +=nuse[i][iboot];
}
