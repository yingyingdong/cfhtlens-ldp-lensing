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
#include "thread_assign.h"

double **gain_xgrid(int NDIV[3],double xrange[2][2]);
void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2]);
void write_grid_num(SHEAR ***shear,int ngroup);
void gain_center(double xc[2],double xrange[2][2]);
void gain_xrange_here(int field,double xrange[2][2]);
void get_line(char *file,int *nobj);
LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
        double xc[2],double mag_limit);
double **D_malloc2_here(int n,int NDIV[2]);
void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xc[2],double mag_limit);
void gain_gstep(double rr,double rmin,double gstep[2],int ntr[2],int flag_lg);
void save_found_zbin(int ****Nfound,int *flag_group,int ngroup,int ntr[2],char *name);
void free_4dhere(int ****Nfound,int ngroup,int nzbin,int ntr0,int *flag_group);
int my_void_here(int nfield,double (*pos)[5],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);
static double step[3];



int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char dir[1024]="/home/fydong/data_dir/data/CFHTLens/";
    char field[4][1024]={"CFHTLens_2015-07-15T01:50:34.tsv","CFHTLens_2015-07-22T03_52_26.tsv",
        "CFHTLens_2015-07-14T02:12:29.tsv","CFHTLens_2015-12-18T03_22_58.tsv"};
	char tmp1[512],tmp2[512],tmp3[512],buf[1024],name[512],file[1024],name_lg[512],
         name_d[512],name_randp[1024]="_";
	int i,j,l,nfit,flagd=1,nfield=0,mask_method=1,load_mask=1,flag_lg=0,flag_randp=1;
    double zbin[6],zmin,zmax,zmeanvoid=0.512,zstep;
    double mag_limit=-20.5,xc[2],mask_ratio=1./10; 
    int nzbin=1,nzbinvoid=1,nmag;	

	//------load data------//
	int nobj = 0,num=13;
	int NDIV[2]={1200,1200},*grid_flag;
    double **arr,xmax[2],ymax[2],xrange[2][2];
	double *zbinlg,*DAzlg,zsteplg; int nzbinlg=1000;
    LENS *Lens;
    nfield = atoi(argv[1]);
    gain_xrange_here(nfield,xrange);
    gain_center(xc,xrange);

    sprintf(file,"%s%s",dir,field[nfield]);
	get_line(file,&nobj);
	arr = D_malloc2(nobj,num);	
	Lens = load_data(file,nobj,num,arr,&nfit,&zmax,xc,mag_limit);	
	free_D_malloc2(arr,nobj);
    zmax=7.;zmin=1e-3;
    zbinlg = zbin_lg(nzbinlg,zmin,zmax,&zsteplg);
    DAzlg  = DA_zlg(nzbinlg,zbinlg);
    printf("zmax = %f \n",zmax);



	//------make grid------//
	double **xgrid;
    int ngroup,ngroup_tmp; double (*group_pos)[5];

	xgrid = D_malloc2_here(2,NDIV);
    step[0] = step[1] = 0.006139;
    gain_NDIV_here(xrange,step,NDIV);
    xgrid = gain_xgrid(NDIV,xrange);
	printf("xmax:%f %f ,ymax:%f %f ,zmax:%f,ndiv0=%d,ndiv1=%d,step=%lg\n",
    xrange[0][0],xrange[0][1],xrange[1][0],xrange[1][1],zmax,NDIV[0],NDIV[1],step[0]);
    
    group_pos = mymalloc(sizeof(double)*5*(2000*2000)*5);

	//------make linklist and make mask------//	
	int *hoc,*ll,**flag_mask;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*nfit);
	
    makell_sub(Lens,nfit,hoc,ll,NDIV,xrange,step,xgrid);
    flag_mask = load_flag_mask(nfield,NDIV[0],NDIV[1]);
		
	//search grid
	int pid,pindex,Nbound=0,kk,zz,nthread=40,*ngroup_assign;
    int ntr[2]={1,60},nuse[nthread],flag_use,ngrp_use=40000;
	double rmax,rmin,rtmp,gstep[2],rcut[5]={1.,1.5,2.,2.5,3.},mul;// /4./60; //10./60
    if(flagd < 1) {rmax=20./60;rmin=0; if(flag_lg==1) rmin=0.1/10;}
    else          {rmax=12000; rmin=0; if(flag_lg==1) rmin=100;}

    if(flag_lg==1) {sprintf(name_lg,"lgrmax%lg-rmin%lg-nrbin%d",rmax,rmin,ntr[1]); }
    else sprintf(name_lg,"rmax%lg-nrbin%d",rmax,ntr[1]);
    if(flagd==1) sprintf(name_d,"adr");
    else sprintf(name_d,"r");
    if(flag_randp == 1) sprintf(name_randp,"randp_");
    gain_gstep(rmax,rmin,gstep,ntr,flag_lg);
    for(zz=0;zz<4;zz++)  //dz
    {
        zstep = 0.05*(1+zz);
        //zstep = 0.1*(zz*0.1+1.6);
        //zstep = 0.1*(1+zz);
        for(kk=0;kk<1;kk++)  //rcut[kk]
        {
            double *random;
            int *flag_group,****Nfound;
            if(flag_randp == 1)
            {
                ngroup=ngrp_use = 100000;
                rand_point2d(group_pos,ngroup,xrange[0][0],xrange[0][1],xrange[1][0],xrange[1][1]);
                for(i=0;i<ngroup;i++)
                {
                    group_pos[i][2] = zmeanvoid;
                    group_pos[i][3] = gain_zbin_lg(DAzlg,group_pos[i][2],zmin,zsteplg);
                }
            }
            else
            {
                ngroup = my_void_here(nfield+1,group_pos,rcut[kk],zstep,zmeanvoid,nzbinvoid,mag_limit,
                         mask_ratio,zmin,zsteplg,DAzlg); //diff here,use the moved ra
            }

            ngroup_assign = thread_assign(nthread,ngroup);
            mul=(ngroup*1./ngrp_use); if(mul<1) mul=1; printf("ngroup=%d,mul=%lg",ngroup,mul); fflush(stdout);
            random = return_drand48(ngroup);
            flag_group = mycalloc(sizeof(int),ngroup); 

            Nfound = mymalloc(sizeof(int ***)*ngroup);
            //for(i=0;i<ngroup;i++)
            //    zero_nfind_thread(Nfound[i],nzbinvoid,ntr);
            //memset(nuse,0,sizeof(int)*nthread);
            
            omp_set_num_threads(nthread); 
            //#pragma omp parallel for schedule(dynamic,1)
            #pragma omp parallel for private(i,rtmp,flag_use)
            for(i=0;i<ngroup;i++)
            {
                if(random[i]<=1./mul) 
                {
                    if(i %10000 == 0) printf("i = %d \n",i);
                    if(flagd==0) rtmp = rmax;
                    else rtmp = (rmax/group_pos[i][3])/au;  //err here
                    if(i<20) printf("rtmp=%lg ",rtmp);
                    if(edge_check(group_pos[i][0],xrange[0][0],xrange[0][1],group_pos[i][1],xrange[1][0],xrange[1][1],rtmp)==1)
                    {
                        Nfound[i] = linklist_search_sub_pos2_countzbin12_record(Lens,group_pos[i],rtmp,rmin,
                               gstep[1],hoc,ll,xrange,step,NDIV,nzbinvoid,zstep,zmeanvoid,flag_mask,
                               flagd,flag_lg,mask_ratio,&flag_use,ntr);
                        if(flag_use == 1) flag_group[i]=1; 
                        else flag_group[i] = 10;
                    }
                }
                //else Nfound[i] = I_malloc3(nzbinvoid,ntr[0],0);
            }
            printf("end \n");


            int tmp=nfield+1;
            sprintf(file,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod%d_step%lg_voidflag.dat",tmp,tmp,mask_method,step[0]);
            //if(load_mask == 0) save_int2d(flag_mask,NDIV[0],NDIV[1],file);
            sprintf(file,"data-2017/number/w%d/w%d%smaskmethod%d_nzbin%d-step%lg-zstep%lg-zmean%lg-rcut%lg_magl%lg-maskratio%lg-%s-%s-mul%d-nfound_ecut.dat",
            //sprintf(file,"data-2017/number/w%d/w%d%smaskmethod%d_nzbin%d-step%lg-zstep%lg-zmin%lg-rcut%lg_magl%lg-maskratio%lg-%s-%s-mul%d-nfound_ecut.dat",
                    tmp,tmp,name_randp,mask_method,nzbin,step[0],zstep,zmeanvoid,rcut[kk],mag_limit,mask_ratio,name_lg,name_d,ngrp_use);
                    //tmp,tmp,name_randp,mask_method,nzbin,step[0],zstep,zmeanvoid-0.5*zstep,rcut[kk],mag_limit,mask_ratio,name_lg,name_d,ngrp_use);
            save_found_zbin(Nfound,flag_group,ngroup,ntr,file);
            fflush(stdout); free(random);
            free_4dhere(Nfound,ngroup,nzbinvoid,ntr[0],flag_group); free(flag_group); 
        }
    }
    free(Lens);
    return 0;
}

void free_4dhere(int ****Nfound,int ngroup,int nzbin,int ntr0,int *flag_group)
{
    int i,j;
    
    for(i=0;i<ngroup;i++) if(flag_group[i]>0)
    {
        for(j=0;j<nzbin;j++)
        {
            free_I_malloc2(Nfound[i][j],ntr0);
        }
        free(Nfound[i]);
    }
    free(Nfound);
}

void gain_NDIV_here(double xrange[2][2],double step[2],int nx[2])
{
    int i,tmp;
    for(i=0;i<2;i++)
        //nx[i] = ceil((xrange[i][1]-xrange[i][0])/step[i]);
        nx[i] = floor((xrange[i][1]-xrange[i][0])/step[i]);  //wrong,but not change because the mask flag also use this method to produce
}

void gain_gstep(double rr,double rmin,double gstep[2],int ntr[2],int flag_lg)
{
    gstep[0] = 2*pi/(ntr[0]);
    if(flag_lg==0) gstep[1] = rr/ntr[1];
    else gstep[1] = (log10(rr)-log10(rmin))/ntr[1];
}


void save_found_zbin(int ****Nfound,int *flag_group,int ngroup,int ntr[2],char *name)
{
    FILE *fp;
    int i,j,k;

    myfopen(fp,name,"w");
    for(k=0;k<ngroup;k++)
    {
        if(flag_group[k]>0)
            for(i=0;i<ntr[0];i++)
            {
                for(j=0;j<ntr[1];j++)
                    fprintf(fp,"%d ",Nfound[k][0][i][j]);
                fprintf(fp,"%d\n",flag_group[k]);
            }
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

int my_void_here(int nfield,double (*pos)[5],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg)
{
    int nline,i,nz2; double tmp;
    FILE *fp;
    char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/w",file[1024];

    sprintf(file,"%s%d/w%d_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmean%lg-rmax%lg_magl%lg-maskratio%lg-pos.dat",
    //sprintf(file,"%s/4_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmin%lg-rmax%lg_magl%lg-maskratio%lg-pos-randompos_mul1.dat",
            dir,nfield,nfield,nzbin_v,zstep,zmean_v,rcut/60.,mag_limit,mask_ratio);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg %lg",&tmp,&pos[i][0],&pos[i][1],&pos[i][2]); //for void(different here from void lensing!!!)
        //fscanf(fp,"%lg %lg %lg",&pos[i][0],&pos[i][1],&pos[i][2]);        //random for void
        //pos[i][2] = zmean_v-0.5*zstep;                                                 
        nz2 = floor( (log10(1+pos[i][2])-log10(1+zmin))/zsteplg );
        pos[i][3] = DAzlg[nz2];
    }
    fclose(fp);
    return nline;
}



