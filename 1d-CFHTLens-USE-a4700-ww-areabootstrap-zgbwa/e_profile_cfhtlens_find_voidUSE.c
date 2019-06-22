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
#include "mask_region.h"
#include "output.h"

double **gain_xgrid(int NDIV[3],double xrange[2][2]);
void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2]);
void write_grid_num(SHEAR ***shear,int ngroup);
void gain_center(double xc[2],double xrange[2][2]);
void gain_xrange_here(int field,double xrange[2][2]);
void gain_group_pos(int NDIV[3],double **xgrid,double (*group_pos)[4]);
void get_line(char *file,int *nobj);
LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
        double xc[2],double mag_limit);
double **D_malloc2_here(int n,int NDIV[2]);
void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xc[2],double mag_limit);
static double step[3];



int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char dir[1024]="/home/dfy/fy/sjtu/data/CFHTLens/";
    char field[4][1024]={"CFHTLens_2015-07-15T01:50:34.tsv","CFHTLens_2015-07-22T03_52_26.tsv",
        "CFHTLens_2015-07-14T02:12:29.tsv","CFHTLens_2015-12-18T03_22_58.tsv"};
	char tmp1[512],tmp2[512],tmp3[512],buf[1024],name[512],file[1024];
    int i,j,nfit,flagd=0,nfield=0,mask_method=1,load_mask=1,nzbin=1;
    double zbin[6],zmean=0.512,zstep,mag_limit=-21.0,xc[2],mask_ratio=1./10; //zmean=0.512	
	
	//------load data------//
	int nobj = 0,num=13;
	int NDIV[2]={1200,1200},*grid_flag;
    double **arr,zmax,xmax[2],ymax[2],xrange[2][2];
	LENS *Lens;
    nfield = atoi(argv[1]); printf("nfield=%d\n",nfield); fflush(stdout);
    gain_xrange_here(nfield,xrange);
    gain_center(xc,xrange);

    sprintf(file,"%s%s",dir,field[nfield]);
	get_line(file,&nobj);
	//nobj = 22724706;
	arr = D_malloc2(nobj,num);	
	Lens = load_data(file,nobj,num,arr,&nfit,&zmax,xc,mag_limit);	
	free_D_malloc2(arr,nobj);
	printf("zmax = %f \n",zmax);



	//------make grid------//
	double **xgrid;
    int ngroup,ngroup_tmp; double (*group_pos)[4];

	xgrid = D_malloc2_here(2,NDIV);
    step[0] = step[1] = 0.006139;
    gain_NDIV_here(xrange,step,NDIV);
    xgrid = gain_xgrid(NDIV,xrange);
	printf("xmax:%f %f ,ymax:%f %f ,zmax:%f,ndiv0=%d,ndiv1=%d,step=%lg\n",
    xrange[0][0],xrange[0][1],xrange[1][0],xrange[1][1],zmax,NDIV[0],NDIV[1],step[0]);
    
    ngroup = NDIV[0]*NDIV[1];
    group_pos = mymalloc(sizeof(double)*4*ngroup);
    gain_group_pos(NDIV,xgrid,group_pos);


	//------make linklist and make mask------//	
	int *hoc,*ll,**flag_mask;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*nfit);
	
    //for(i=0;i<2;i++) {step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
    //    printf("step[%d]=%f\n",i,step[i]);}

    makell_sub(Lens,nfit,hoc,ll,NDIV,xrange,step,xgrid);
    //if(load_mask == 0) flag_mask = mask_region(NDIV[0],NDIV[1],hoc,ll,step,xgrid,xc[0],xc[1],mask_method,nfield);
    flag_mask = load_flag_mask(nfield,NDIV[0],NDIV[1]);
		
	//search grid
	int **Nfound,pid,pindex,Nbound=0,kk,zz;
	double rmax[5]={1./60,1.5/60,2./60,2.5/60,3./60},rtmp;// /4./60; //10./60
	double rmin = 0.;//0.1/60;
    Nfound   = mymalloc(sizeof(int *)*ngroup);


    for(zz=2;zz<4;zz++)
    {
        zstep = 0.05*(1+zz);
        //zstep = 0.1*(zz*0.1+1.6);
        for(kk=0;kk<5;kk++)
        {
            omp_set_num_threads(80);
            //#pragma omp parallel for schedule(dynamic,1)
#pragma omp parallel for private(i,rtmp)
            for(i=0;i<ngroup;i++)
            {
                if(i %10000 == 0) printf("i = %d \n",i);
                if(flagd==0) rtmp = rmax[kk];
                else rtmp = (rmax[kk]/group_pos[i][3])/au;  //err here
                Nfound[i] = linklist_search_sub_pos2_countzbin12(Lens,group_pos[i],rtmp,
                        hoc,ll,xrange,step,NDIV,nzbin,zstep,zmean,flag_mask,flagd,mask_ratio);
            }
            printf("end \n");



            int tmp=nfield+1;
            sprintf(file,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_step%lg_voidflag.dat",tmp,tmp,mask_method,step[0]);
            //if(load_mask == 0) save_int2d(flag_mask,NDIV[0],NDIV[1],file);
            sprintf(file,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_nzbin%d-step%lg-zstep%lg-zmean%lg-rmax%lg_magl%lg-maskratio%lg-nfound.dat",tmp,tmp,mask_method,nzbin,step[0],zstep,zmean,rmax[kk],mag_limit,mask_ratio);
            sprintf(tmp1,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_nzbin%d-step%lg-zstep%lg-zmean%lg-rmax%lg_magl%lg-maskratio%lg-pos.dat",tmp,tmp,mask_method,nzbin,step[0],zstep,zmean,rmax[kk],mag_limit,mask_ratio);
            //sprintf(file,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_nzbin%d-step%lg-zstep%lg-zmin%lg-rmax%lg_magl%lg-found.dat",tmp,tmp,mask_method,nzbin,step[0],zstep,zmin,rmax[kk],mag_limit);
            //sprintf(tmp1,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_nzbin%d-step%lg-zstep%lg-zmin%lg-rmax%lg_magl%lg-pos.dat",tmp,tmp,mask_method,nzbin,step[0],zstep,zmin,rmax[kk],mag_limit);
            save_found_zbin_o(Nfound,nzbin,NDIV[0],NDIV[1],file);
            save_void_pos_o(Nfound,flag_mask,NDIV,nzbin,group_pos,tmp1,zstep,zmean,xc);
        }
    }
    free(Lens);
    return 0;
}

void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2])
{
    int i,tmp;
    for(i=0;i<2;i++)
        //nx[i] = ceil((xrange[i][1]-xrange[i][0])/step[i]);
        nx[i] = floor((xrange[i][1]-xrange[i][0])/step[i]);  //wrong
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

void write_grid_num(SHEAR ***shear,int ngroup)
{
    int i,nx;
    FILE *fp;

    nx = (int)(sqrt(ngroup)); printf("nx=%d\n",nx);
    myfopen(fp,"num.dat","w");
    for(i=0;i<ngroup;i++)
    {
        fprintf(fp,"%d ",shear[i][0][0].nrr);
        if((i+1)%nx==0) fprintf(fp,"\n");
    }
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


LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
		double xc[2],double mag_limit)
{
	int i,n,tmp = 0;
	FILE *fp;          
	char chararr[512];
	LENS *Lens;
	myfopen(fp,file,"r");                                                               
	fscanf(fp,"%*[^\n]%*c");  //skip the first line                                                          
	n = 0;                                                                            
	while(n < nobj)                                                                
	//while(n < 5)                                                                
	{                                                                                                                                            
		for(i=0;i<3;i++) fscanf(fp,"%s",chararr);	
		for(i=0;i<num;i++)                                                               
			fscanf(fp,"%lf",&arr[n][i]);                                                    
		//if(arr[n][7] == 0 && arr[n][8] <= 1.2 && arr[n][12] > -99) tmp++;
        if(arr[n][11] > -98 && arr[n][11]<98 && arr[n][11]<mag_limit)
            tmp ++;
		n ++;		                                                                  
	}                                                                  
	fclose(fp);                    
	*nfit = tmp;

	printf("number of good data is %d \n",tmp);
	Lens = mymalloc(sizeof(LENS)*tmp);
	good_data(arr,nobj,tmp,Lens,zmax,xc,mag_limit);
	return Lens;
}

void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xc[2],double mag_limit)
{
	int i,tmp = 0,flag;
	float ymax = -10,xmax = 0,zmax = 0,xmin = 1000,ymin = 1000;
	for(i=0;i<nobj;i++)
		//if(arr[i][7] == 0 && arr[i][8] <= 1.2 && arr[i][12] > -99)
		if(arr[i][11] > -98 && arr[i][11]<98 && arr[i][11]<mag_limit)
        {
			Lens[tmp].e1 = arr[i][4];  
			Lens[tmp].e2 = arr[i][5]-arr[i][10];
			Lens[tmp].m  = arr[i][9];
			Lens[tmp].mag = arr[i][11];  //absolute L
			//Lens[tmp].mag = arr[i][12];  //mag
			Lens[tmp].wt = arr[i][6]*0.01;
			Lens[tmp].pos[0] = (-arr[i][0]-xc[0])*cos(au*arr[i][1])+xc[0];
			Lens[tmp].pos[1] = arr[i][1];
			Lens[tmp].pos[2] = arr[i][8];
			
			if(xmax < arr[i][0]) xmax = arr[i][0];
			if(ymax < arr[i][1]) ymax = arr[i][1];
			if(zmax < arr[i][8]) zmax = arr[i][8];
			if(xmin > arr[i][0]) xmin = arr[i][0];
			if(ymin > arr[i][1]) ymin = arr[i][1];
						
			tmp ++;
		}
	printf("zmax = %f \n",zmax);
	*z = zmax;
}




