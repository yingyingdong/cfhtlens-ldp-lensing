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

void write_grid_num(SHEAR ***shear,int ngroup);
void gain_center(double xc[2],double xrange[2][2]);
void gain_xrange_here(int field,double xrange[2][2]);
void gain_group_pos(int NDIV[3],double **xgrid,double (*group_pos)[4]);
void gain_random_for_void(int **flag_mask,double (*group_pos)[5],int ngroup,double xrange[2][2],
        int NDIV[2],double step[2],char *dir);
void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2]);
static double step[3];
int my_void_here(int field,double (*pos)[5],double rcut,double zstep,double zmin_v,int nzbin_v,
                double mag_limit,double mask_ratio,char *file);
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
	char tmp1[512],dir[512],name[512],file[1024],file_flag[1024];
	int i,j,nfit,flagd=0,nfield=3,mask_method=1;
	
	//------load data------//
	int NDIV[2]={1200,1200},*grid_flag;
	double xrange[2][2],**xgrid,xc[2],rmax = 1./60;

    gain_xrange_here(nfield,xrange);
	printf("xmax:%f %f ,ymax:%f %f\n",xrange[0][0],xrange[0][1],
			xrange[1][0],xrange[1][1]);
    //for(i=0;i<2;i++) {step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
    //    printf("step[%d]=%f\n",i,step[i]);}
    step[0] = step[1] = 0.006139;
    gain_NDIV_here(xrange,step,NDIV);

    int **flag_mask,tmp,ngroup;
    double (*group_pos)[5],rcut[5]={1.,1.5,2.,2.5,3.};
    double zstep=0.1,zminvoid=0.4,mag_limit=-20.5,mask_ratio=0.1; int kk,nzbinvoid=1;
    tmp = nfield+1;
    group_pos = mymalloc(sizeof(double)*5*(1200*1200));
    for(kk=0;kk<5;kk++)
    {
        ngroup = my_void_here(tmp,group_pos,rcut[kk],zstep,zminvoid,nzbinvoid,mag_limit,mask_ratio,file);

        //sprintf(dir,"data-2017/mask_flag/maskmethod%d_nlimit12_%d_rmax%lg/w%d",mask_method,NDIV[0],rmax,tmp);
        sprintf(dir,"data-2017/mask_flag/maskmethod1_USE/w%d",tmp);
        //sprintf(file_flag,"%s/w%d_flag_maskmethod%d_ngrid%d_void.dat",dir,tmp,mask_method,NDIV[0]); //flag
        sprintf(file_flag,"%s/w%d_maskmethod1_step0.006139_voidflag.dat",dir,tmp); //flag
        flag_mask = load_int2d(NDIV[0],NDIV[1],file_flag);
        gain_random_for_void(flag_mask,group_pos,ngroup,xrange,NDIV,step,file);
    }
    return 0;
}


void gain_random_for_void(int **flag_mask,double (*group_pos)[5],int ngroup,double xrange[2][2],
        int NDIV[2],double step[2],char *file)
{
    int i,j,k,nx[2],mul=1,num;
    double x[2],(*xrand)[4];
    char ofile[1024];
    xrand = mymalloc(sizeof(double)*4*ngroup*mul);

    srand((unsigned)time(NULL));
    num = 0;
    for(i=0;i<ngroup;i++)
    {
        for(k=0;k<mul;k++)
        {
            do{
                for(j=0;j<2;j++)
                {
                    x[j]  = (rand()/((double)RAND_MAX+1.))*(xrange[j][1]-xrange[j][0])+xrange[j][0];
                    nx[j] = floor((x[j]-xrange[j][0])/step[j]);
                }
            }while(nx[0]>=NDIV[0] && nx[1]>=NDIV[1] && flag_mask[nx[1]][nx[0]]==0);
            for(j=0;j<2;j++) xrand[num][j]=x[j];
            xrand[num][2] = group_pos[i][2];
            num++;
        }
    }
    
    //sprintf(file,"%s/random_pos_mul%d.dat",dir,mul);
    sprintf(ofile,"%s-randompos_mul%d.dat",file,mul); 
    save_double2d4(xrand,3,ngroup*mul,ofile);
}

int my_void_here(int field,double (*pos)[5],double rcut,double zstep,double zmin_v,int nzbin_v,
        double mag_limit,double mask_ratio,char *file)
{
    int nline,i,nz2; double tmp;
    FILE *fp;
    char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/",infile[1024];

    sprintf(file,"%sw%d/w%d_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmin%lg-rmax%lg_magl%lg-maskratio%lg-pos",
            dir,field,field,nzbin_v,zstep,zmin_v,rcut/60.,mag_limit,mask_ratio);
    sprintf(infile,"%s.dat",file);
    printf("%s\n",infile);
    GET_line(infile,&nline);
    myfopen(fp,infile,"r");
    for(i=0;i<nline;i++)
        fscanf(fp,"%lg %lg %lg %lg",&pos[i][0],&tmp,&pos[i][1],&pos[i][2]);
    fclose(fp);
    return nline;
}


void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2])
{
    int i,tmp;
    for(i=0;i<2;i++)
        //nx[i] = ceil((xrange[i][1]-xrange[i][0])/step[i]);
        nx[i] = floor((xrange[i][1]-xrange[i][0])/step[i]); //wrong
    printf("ndiv:%d,%d\n",nx[0],nx[1]);
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



#undef omegal0 
#undef omegam0 
#undef H0 
#undef DH 
#undef au 
