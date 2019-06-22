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
#pragma comment(lib, "libgsl_d.lib")  
#pragma comment(lib, "libgslcblas_d.lib") 
#include "func.h"
#include "struct.h"
#include "Romberg.h"
#include "power.h"
#include "clusters.h"
#include "get_z_dal.h"
#include "thread_assign.h"
#include "output.h"

void get_line(char *file,int *nobj);
double **D_malloc2_here(int n,int NDIV[2]);
float arcsin(double x,double y);
static double step[2];
void check1d(int *arr,int n);
void check_shear(SHEAR **shear,int n1,int n2);

void free_arr_use(int **found,int *Nfound,int ngroup,double *random);
void zero_shear_here(SHEAR ***shear,SHEARR **shearr,int nzbin,int ntr[2]);
//double angu_dia_d(double z)
//{
//	return 1./sqrt(omegal0+omegam0*(1+z)*(1+z)*(1+z));
//}


int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char file[1024]="/home/dfy/data/CFHTLens/4wpz.tsv";
    //char file[1024]="/home/fydong/data_dir/data/CFHTLens/CFHTLens_2015-07-15T01:50:34-pz-jz-fields.tsv";
	char tmp1[512],tmp2[512],tmp3[512],wtname[512],dirout[1024];
	char buf[1024],name[512],tmpwt[512]; 
    double magl[4]={-20.5,-21,-21.5,-22.},mag_limit,zmeanvoid=0.35; //zmeanvoid=0.512,0.45
	int i,j,l,nfit,flagc=-1,flagd=1,flag_lg=1,wt_sigma=1,flag_com=1,flag_xy=0,flag_density=1,nmag;
    if(flag_density == 0) sprintf(dirout,"onetime-pz");
    if(flag_density == 1) sprintf(dirout,"onetime-sigma-pz");
    //mag_limit = atoi(argv[1]);
	
	//------load data------//
	int nobj = 0,num=15,nzbinlg=1000;
	double **arr,*zbinlg,*DAzlg,zsteplg;
	double zmax=7.,zmin=1e-3,xmax[2],ymax[2];
	LENS *Lens;

    //printf("%f\n",sigma(0.39,0.65,1e-3,zsteplg,0,DAzlg));
	get_line(file,&nobj);
	//nobj = 22724706;
	arr = D_malloc2(nobj,num);	
	Lens = load_data_o(file,nobj,num,arr,&nfit,&zmax,xmax,ymax,zmeanvoid-0.1);	
	free_D_malloc2(arr,nobj);
    zbinlg = zbin_lg(nzbinlg,zmin,zmax,&zsteplg);
    DAzlg  = DA_zlg_ww(nzbinlg,zbinlg);
    printf("data is load,zmax = %f \n",zmax); fflush(stdout);



	//------make grid------//
	int NDIV[2]={4096,4096};
    double xrange[2][2];
	double **xgrid;
	xgrid = D_malloc2_here(2,NDIV);
	for(i=0;i<2;i++)
		xrange[1][i] = ymax[i];
	xrange[0][0] = -xmax[1];  
	xrange[0][1] = -xmax[0];
	
	printf("xmax:%f %f ,ymax:%f %f ,zmax:%f\n",xrange[0][0],xrange[0][1],
			xrange[1][0],xrange[1][1],zmax);  fflush(stdout);

	for(i=0;i<2;i++) 
		for(j=0;j<=NDIV[i];j++)
			xgrid[i][j] = xrange[i][0]+j*(xrange[i][1]-xrange[i][0])/NDIV[i];

	
	//------make linklist------//	
	int *hoc,*ll;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*nfit);
	
	for(i=0;i<2;i++) step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
	makell_sub(Lens,nfit,hoc,ll,NDIV,xrange,step,xgrid); 
	printf("ll is done\n");fflush(stdout);
	//check_num(Lens,NDIV,xgrid,hoc,ll);
		
	//search grid
	//first search source in zgrid
	int zgrid1 = 0,zgrid2 = zgrid1+1,dz=3;  //zgird = 2
	double zbin_c[2]={0.2,0.58},pz=0.2;
	//double zbin_b[5][2]={{0.58,0.72},{0.72,0.86},{0.86,1.02},{1.02,1.3},{1.3,2.5}};
	double zbin_b[1][2]={{0.58,2.5}};
	int nzbin=1,ngroup;
    //double zbin[2]={0.1,0.15};
	double (*group_pos)[6];
    group_pos = mymalloc(sizeof(double)*6*(2000*2000));

	int pid,pindex,Nbound=0;
	double mul,rmax,rtmp,rmin = 100.;
    double gstep[2],rcut[6]={1.,1.5,2.,2.5,3.}; //arcmin
    //double gstep[2],rcut[6]={1.,1.2,1.5,2.,2.5,3.}; //arcmin
    int ntr[2] = {1,60}; if(flagd < 1) ntr[1] = 20;
    int nthread=20;

    gain_nr_gstep_o(ntr,gstep,flag_lg,flagd,flag_xy,&rmin,&rmax);
    SHEAR ****shear,***tmp_shear;
    SHEARR ***shearr;
    tmp_shear = mymalloc(sizeof(SHEAR **)*nthread);
    for(i=0;i<nthread;i++) tmp_shear[i] = shear_calloc(ntr[0],ntr[1]);
    shear  = shear_thread_assign_nboot(nthread,ntr);
    shearr = shearr_thread_assign_nboot(nthread,ntr);
        
    double zstep,mask_ratio=0.1; 
    int kk,iboot,nzbinvoid=1,nn,**nuse,zz,ngroup_use,ngrp_use;
    nuse = I_calloc2(nthread,nboot+1);
    for(nmag=1;nmag<4;nmag++) 
    {
        mag_limit = magl[nmag];
        for(zz=3;zz<4;zz++)
        {
            zstep = 0.05*(1+zz);
            for(kk=0;kk<2;kk++)
            {
                ngroup = my_void_4field_groupflag(group_pos,rcut[kk],zstep,zmeanvoid,nzbinvoid,mag_limit,mask_ratio
                        ,zmin,zsteplg,DAzlg);
                int **found,*Nfound,*ngroup_assign; 
                FLAG_GROUP *flag_group;
                ngroup_assign = thread_assign(nthread,ngroup);
                //ngrp_use = MIN_I(ngroup_use,ngroup); 
                ngrp_use = ngroup;
                flag_group = flag_group_malloc_areaboot(nboot,ngroup,group_pos);

                zero_shear_thread_nboot(shear,shearr,nuse,ntr,nthread);
                omp_set_num_threads(nthread);
                //#pragma omp parallel for schedule(dynamic,1)
#pragma omp parallel for private(l,i,rtmp)
                for(l=0;l<nthread;l++)
                {
                    for(i=ngroup_assign[l];i<ngroup_assign[l+1];i++)
                        if(flag_group[nboot].flag[i]>0)
                        {
                            if(flagd==0) rtmp = rmax/cos(au*group_pos[i][1]);
                            else rtmp = (rmax/group_pos[i][3])/au/cos(au*group_pos[i][1]);
                            if(flag_xy == 0) 
                                linklist_search_sub_pos2_pz_thread(Lens,group_pos[i],tmp_shear[l],rtmp,rmin,hoc,ll,xrange,step,gstep,zbin_b[0],ntr,zsteplg,\
                                        NDIV,nfit,DAzlg,zmin,flagd,flag_lg,wt_sigma,flag_com,flag_density,pz);
                            add_shear_to_iboot(tmp_shear[l],shear[l],flag_group,nuse[l],ntr,i);
                        }
                }
                printf("end \n");


                //for(i=0;i<ngroup;i++) {if(i%10==0)nn=0;else nn=10; found[i] = mymalloc(sizeof(int)*nn);}
                gain_file_name_xy(flag_xy,flag_lg,flagd,flagc,wt_sigma,flag_com,tmp1,tmp2,tmpwt,ngrp_use,
                        mag_limit,zmeanvoid,rcut[kk],nzbinvoid,zstep);

                for(iboot=0;iboot<=nboot;iboot++)
                {
                    add_shear_thread(nthread,shear,ntr,iboot);
                    add_nuse_thread(nuse,nthread,iboot);
                    sprintf(file,"data-2017/%s/4jzfield-ww%lg-%s-%d-%d-shear-cali2-%scra-%lg-%lg-%s-pz%lgomp-iboot%d-group%d-%d-%d.dat",
                            dirout,ww,tmp2,ntr[0],ntr[1],tmp1,zbin_b[0][0],zbin_b[0][1],tmpwt,pz,iboot,3000,stepdec,stepra); 
                    printf("filename:%s\n",file); fflush(stdout);
                    if(flag_xy == 0)
                        shear_output_r(ntr,shear[0][iboot],shearr[0][iboot],file,gstep,rmin,flag_lg,nuse[0][iboot]);
                    else
                        shear_output_xy(ntr,shear[0][iboot],shearr[0][iboot],file,gstep,rmax/sqrt(2.));	
                }
                free(ngroup_assign);
                for(i=0;i<=nboot;i++) free(flag_group[i].flag); free(flag_group);
            }
        }
    }
	return 0;
}




void zero_shear_here(SHEAR ***shear,SHEARR **shearr,int nzbin,int ntr[2])
{
    int i,j;
    for(i=0;i<nzbin;i++)
    {
        for(j=0;j<ntr[0];j++)
            memset(shear[i][j],0,sizeof(SHEAR)*ntr[1]);
        memset(shearr[i],0,sizeof(SHEARR)*ntr[1]);
    }
}





void free_arr_use(int **found,int *Nfound,int ngroup,double *random)
{
    int i;
    free(random);
    free_I_malloc2(found,ngroup);
    free(Nfound);
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

	//fscanf(fp,"%*[^\n]%*c");
	while(fgets(buf,1024,fp) != NULL)
	{
		if(buf[strlen(buf) - 1] == '\n')
			n++;
	}
	fclose(fp);
	*nobj = n;
	printf("nobj = %d \n",*nobj);

}




int if_range(double xrange[2][2],double tmp[2])
{
	if(-tmp[0]>xrange[0][0] && -tmp[0]<xrange[0][1]
			&& tmp[1]>xrange[1][0] && tmp[1]<xrange[1][1])
		return 1;
	else return 0;
	
}



void check1d(int *arr,int n)
{
	int i;
	for(i=0;i<n;i++)
		if(arr[i] != 0) printf("use of memset is wrong!\n");
}




