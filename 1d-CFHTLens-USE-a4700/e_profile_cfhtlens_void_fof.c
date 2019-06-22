#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fftw3.h"
#include "omp.h"

#include "func.h"
#include "struct.h"
#include "Romberg.h"
#include "power.h"
#include "clusters.h"
#include "get_z_dal.h"
#include "fft_m.h"

double **gain_xgrid(int NDIV[3],double xrange[2][2]);
void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2]);
void gain_xrange_here(int field,double xrange[2][2]);
void gain_group_pos(int NDIV[3],double **xgrid,double (*group_pos)[4]);
void get_line(char *file,int *nobj);
double **D_malloc2_here(int n,int NDIV[2]);
static double step[3];

LENS *copy_to_Lens(int ngroup,double (*group_pos)[4]);
int judge_fof(LENS *Lens,int nfound,int *found,int flag[30000]);
void flag_fof(int ntmp,LENS *Lens,int nfound,int *found);
void multi_fof(int flag_tmp,int flag[30000],LENS *Lens,int Nfound,int *found,int ngroup);
void mean_fof_pos(LENS *Lens,double pos[][3],int ntmp,int ngroup,double rmax,int nfield);
void check_num_here(int NDIV[2],int *hoc,int *ll);

#define omegal0 0.7
#define omegam0 0.3
#define H0 72.
#define au (pi/180.)  //error
//double angu_dia_d(double z)
//{
//	return 1./sqrt(omegal0+omegam0*(1+z)*(1+z)*(1+z));
//}


int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char tmp1[512],tmp2[512],tmp3[512],buf[1024],name[512],file[1024];
	int i,j,flagc=-1,flagd=0,nfield=3;
	int NDIV[2]={100,100},*grid_flag;
    double xrange[2][2],**xgrid; LENS *Lens;
    gain_xrange_here(nfield,xrange);

	//------make grid------//
    int ngroup,ngroup_tmp; double (*group_pos)[4];

	xgrid = D_malloc2_here(2,NDIV);
    for(i=0;i<2;i++) step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
	printf("xmax:%f %f ,ymax:%f %f ,ndiv0=%d,ndiv1=%d,step=%lg\n",
    xrange[0][0],xrange[0][1],xrange[1][0],xrange[1][1],NDIV[0],NDIV[1],step[0]);
    
    group_pos = mymalloc(sizeof(double)*4*30000);
    ngroup = read_cluster(flagc,group_pos); printf("ngroup=%d\n",ngroup);
    Lens = copy_to_Lens(ngroup,group_pos);

	//------make linklist and make mask------//	
	int *hoc,*ll,**flag_mask;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*ngroup);
    makell_sub(Lens,ngroup,hoc,ll,NDIV,xrange,step,xgrid);	
	check_num_here(NDIV,hoc,ll);
    //search grid
	int **found,*Nfound,**fof,*nfof,pid,pindex,Nbound=0,ntmp,flag[30000]={0},flagh,n3=30000;
	double rmax,rtmp,rmin=0.;// /4./60; //10./60
    rmax = 1./60;//MAX_D(step[0],step[1]);//5./60;
    found = mymalloc(sizeof(int *)*ngroup);
    Nfound   = mycalloc(sizeof(int),ngroup);
    fof  = I_calloc2(ngroup,200);
    nfof = mycalloc(sizeof(int),ngroup);

    ntmp = 1;
	for(i=0;i<ngroup;i++) 
	{
		if(flagd==0) rtmp = rmax;
		else rtmp = (rmax/group_pos[i][3])/au;
        found[i] = linklist_search_sub_pos2(Lens,group_pos[i],rtmp,
                &Nfound[i],hoc,ll,xrange,step,NDIV,ngroup);
        memset(flag,0,sizeof(int)*n3);
        flagh = judge_fof(Lens,Nfound[i],found[i],flag);
        printf("nfound[%d]=%d,flah=%d\n",i,Nfound[i],flagh);
        if(flagh == 0) 
        {
            flag_fof(ntmp,Lens,Nfound[i],found[i]);
            ntmp ++;
        }
        else if(flagh == 1) //has some pos already in some fof
            flag_fof(flag[0],Lens,Nfound[i],found[i]);
        else
            multi_fof(flagh,flag,Lens,Nfound[i],found[i],ngroup);
    }
	printf("end \n");

    double pos[ntmp][3];
    mean_fof_pos(Lens,pos,ntmp,ngroup,rmax,nfield+1);

	
	return 0;
}

void check_num_here(int NDIV[2],int *hoc,int *ll)
{
    int i,j,index,pid,tmp=0;

    for(i=0;i<NDIV[0];i++)
        for(j=0;j<NDIV[1];j++)
        {
            index = i+j*NDIV[0];
            pid = hoc[index];
            while(pid >= 0)
            {
                pid = ll[pid];
                tmp ++;
            }
        }
    printf("check num=%d\n",tmp);
}

int judge_fof(LENS *Lens,int nfound,int *found,int flag[30000])
{
    int i,j,pid,flagh,fflag;
    
    flagh = 0;
    for(i=0;i<nfound;i++)
    {
        pid = found[i]; 
        fflag = 0;
        if(Lens[pid].flag>0)
        {
            for(j=0;j<flagh;j++)
                if(Lens[pid].flag==flag[j]) {fflag=1;break;}
            if(fflag==0) {flag[flagh]=Lens[pid].flag;flagh ++;}
        }
    }
    return flagh;
}

void flag_fof(int ntmp,LENS *Lens,int nfound,int *found)
{
    int i,pid;

    for(i=0;i<nfound;i++)
    {
        pid = found[i]; 
        //printf("%d ",pid);
        Lens[pid].flag = ntmp;
    }
    //printf("\n");
}

void multi_fof(int flagh,int flag[30000],LENS *Lens,int nfound,int *found,int ngroup)
{
    int i,j;

    flag_fof(flag[flagh-1],Lens,nfound,found);
    
    for(j=0;j<ngroup;j++)
    {
        for(i=0;i<flagh-1;i++)
            if(Lens[j].flag==flag[i]) Lens[j].flag = flag[flagh-1];
    }

}

void mean_fof_pos(LENS *Lens,double pos[][3],int ntmp,int ngroup,double rmax,int nfield)
{
    int i,j,k,num[ntmp],temp=0,nmin=4;
    double ***POS;
    FILE *fp1,*fp2; char dir[512],file1[1024],file2[1024];

    sprintf(dir,"data-2017/mask_flag/maskmethod1_USE/w%d",nfield);
    for(i=0;i<ntmp;i++) for(j=0;j<3;j++) pos[i][j]=0;
    POS = mymalloc(sizeof(double **)*ntmp); for(i=0;i<ntmp;i++) POS[i] = D_calloc2(2000,3);
    memset(num,0,sizeof(int)*ntmp);
    for(i=0;i<ngroup;i++)
    {
        for(k=1;k<=ntmp;k++) if(Lens[i].flag == k)
        {
            for(j=0;j<3;j++)
            {
                pos[k-1][j] +=Lens[i].pos[j];
                POS[k-1][num[k-1]][j] = Lens[i].pos[j];
            }
            num[k-1] ++;
            temp ++;
        }
    }
    printf("ntmp=%d\n",ntmp);
    printf("total num=%d\n",temp);
    sprintf(file1,"%s/w%d_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.05_magl-20-fofr%lgpos.dat",dir,nfield,rmax*60);
    sprintf(file2,"%s/w%d_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.05_magl-20-fofr%lgPOS-n%d.dat",dir,nfield,rmax*60,nmin);
    myfopen(fp1,file1,"w");
    myfopen(fp2,file2,"w");
    for(k=0;k<ntmp;k++) if(num[k]>0)   
    {
        for(j=0;j<3;j++)
            fprintf(fp1,"%lf ",pos[k][j]/num[k]);
        fprintf(fp1,"%d\n",num[k]);
        if(num[k] >= nmin) for(i=0;i<num[k];i++)
        {
            for(j=0;j<3;j++)
                fprintf(fp2,"%lf ",POS[k][i][j]);
            fprintf(fp2,"%d\n",k); //the k'th fof
        }
    }
    fclose(fp1);
    fclose(fp2);
}


LENS *copy_to_Lens(int ngroup,double (*group_pos)[4])
{
    int i,j;
    LENS *Lens;
    Lens = mycalloc(sizeof(LENS),ngroup);
    for(i=0;i<ngroup;i++)
    {
        for(j=0;j<3;j++)
        {
            Lens[i].pos[j] = group_pos[i][j];
            printf("%lg ",Lens[i].pos[j]);
        }
        printf("\n");
    }
    return Lens;

}

void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2])
{
    int i,tmp;
    for(i=0;i<2;i++)
        //nx[i] = ceil((xrange[i][1]-xrange[i][0])/step[i]);
        nx[i] = floor((xrange[i][1]-xrange[i][0])/step[i]);  //wrong
}

void save_found_zbin(int **Nfound,int nzbin,int nx,int ny,char *name)
{
    FILE *fp;
    int i,j,k;

    myfopen(fp,name,"w");
    for(k=0;k<nzbin;k++)
        for(i=0;i<nx*ny;i++)
        {
            fprintf(fp,"%d ",Nfound[i][k]);
            if((i+1)%nx==0) fprintf(fp,"\n");
        }
    fclose(fp);
}


double **gain_xgrid(int NDIV[3],double xrange[2][2])
{
    int i,j;  double **xgrid;
    FILE *fp;
    xgrid = D_malloc2_here(2,NDIV);

    for(i=0;i<2;i++)
    {
        for(j=0;j<=NDIV[i];j++)
        {
            xgrid[i][j] = xrange[i][0]+j*(xrange[i][1]-xrange[i][0])/NDIV[i];
        }
    }
    return xgrid;
}


void gain_center(double xc[2],double xrange[2][2])
{
    xc[0] = 0.5*(xrange[0][0]+xrange[0][1]);
    xc[1] = 0.5*(xrange[1][0]+xrange[1][1]);
}


void gain_xrange_here(int field,double xrange[2][2])
{
    int i,j;

    if(field == 0)
    {
        xrange[0][0] = -38.820652-0.1; xrange[0][1] = -30.17882+0.1;
        xrange[1][0] = -11.244478-0.1; xrange[1][1] = -3.6772671+0.1;
    }
    if(field == 1)
    {
        xrange[0][0] = -136.84346-0.1; xrange[0][1] = -132.0583+0.1;
        xrange[1][0] = -5.6952291-0.1; xrange[1][1] = -0.952663+0.1;
    }
    if(field == 2)
    {
        xrange[0][0] = -220.38112-0.1; xrange[0][1] = -208.55965+0.1;
        xrange[1][0] = 51.19738-0.1; xrange[1][1] = 57.80508+0.1;
    }
    if(field == 3)
    {
        xrange[0][0] = -335.70685-0.1; xrange[0][1] = -329.97742+0.1;
        xrange[1][0] = -1.029174-0.1; xrange[1][1] = 4.6142831+0.1;
    }
}

void gain_group_pos(int NDIV[3],double **xgrid,double (*group_pos)[4])
{
    int i,j;
    FILE *fp;

    for(i=0;i<NDIV[0]*NDIV[1];i++)
    {
        group_pos[i][0] = 0.5*(xgrid[0][i%NDIV[0]]+xgrid[0][i%NDIV[0]+1]);
        group_pos[i][1] = 0.5*(xgrid[1][i/NDIV[0]]+xgrid[1][i/NDIV[0]+1]);
    }
    //myfopen(fp,"group_pos.dat","w");
    //for(i=0;i<NDIV[1]*NDIV[0];i++)
    //{
    //    fprintf(fp,"%lg %lg\n",group_pos[i][0],group_pos[i][1]);
    //}
    //fclose(fp);
}


double **D_malloc2_here(int n,int NDIV[2])
{
	int i,j;
	double **arr;

	arr = mymalloc(sizeof(double *)*n);
	for(i=0;i<n;i++)
		arr[i] = mymalloc(sizeof(double)*(NDIV[i]+1));
	return arr;
}

void get_line(char *file,int *nobj)
{
	FILE *fp; int n = 0;
	char buf[1024];
	myfopen(fp,file,"r");

	fscanf(fp,"%*[^\n]%*c");
	while(fgets(buf,1024,fp) != NULL)
	{
		if(buf[strlen(buf) - 1] == '\n')
			n++;
	}
	fclose(fp);
	*nobj = n;
	printf("nobj = %d \n",*nobj);

}

#undef omegal0 
#undef omegam0 
#undef H0 
#undef DH 
#undef au 
