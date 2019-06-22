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

void get_line(char *file,int *nobj);
LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
		        double xmax[2],double ymax[2]);
double **D_malloc2_here(int n,int NDIV[2]);
void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2]);
//int clusters(double (*pos)[4],int num,double xrange[2][2]);
float arcsin(double x,double y);
void mean_shear(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
		double gstep[2],int ntr[2],double zbin[2],double rr,double rmin,int flagd,int flag_lg);
double et1(double e1,double e2,double sin2a,double cos2a);
double et2(double e1,double e2,double sin2a,double cos2a);
static double step[2];
void check1d(int *arr,int n);
void check_shear(SHEAR **shear,int n1,int n2);
void shear_output(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],double rr);
void gain_gstep(double rr,double rmin,double gstep[2],int ntr[2],int flag_lg);
void gain_fileo_name(int flagd,int flagc,char *tmp1,char *tmp2);

#define omegal0 0.7
#define omegam0 0.3
#define H0 72.
//#define DH 3000000 //(kpc/h)
#define au (pi/180.)  //error
#define nzbin 1  //error
//double angu_dia_d(double z)
//{
//	return 1./sqrt(omegal0+omegam0*(1+z)*(1+z)*(1+z));
//}


int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char file[1024]="/home/fydong/data_dir/data/CFHTLens/4-jz-fields.tsv";
	char tmp1[512],tmp2[512],tmp3[512];
	char buf[1024],name[512];
	int i,j,nfit,flagc=-1,flagd=0,flag_lg=0;

	
	
	//------load data------//
	int nobj = 0,num=13;
	double **arr;
	double zmax,xmax[2],ymax[2];
	LENS *Lens;
	get_line(file,&nobj);
	//nobj = 22724706;
	arr = D_malloc2(nobj,num);	
	Lens = load_data(file,nobj,num,arr,&nfit,&zmax,xmax,ymax);	
	free_D_malloc2(arr,nobj);
	printf("zmax = %f \n",zmax);



	//------make grid------//
	int NDIV[2]={50,50};
	double xrange[2][2];
	double **xgrid;
	xgrid = D_malloc2_here(2,NDIV);
	for(i=0;i<2;i++)
		xrange[1][i] = ymax[i];
	xrange[0][0] = -xmax[1];  
	xrange[0][1] = -xmax[0];
	
	printf("xmax:%f %f ,ymax:%f %f ,zmax:%f\n",xrange[0][0],xrange[0][1],
			xrange[1][0],xrange[1][1],zmax);

	for(i=0;i<2;i++) 
		for(j=0;j<=NDIV[i];j++)
			xgrid[i][j] = xrange[i][0]+j*(xrange[i][1]-xrange[i][0])/NDIV[i];

	
	//------make linklist------//	
	int *hoc,*ll;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*nfit);
	
	for(i=0;i<2;i++) step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
	makell_sub(Lens,nfit,hoc,ll,NDIV,xrange,step,xgrid);

	//check_num(Lens,NDIV,xgrid,hoc,ll);
		
	//search grid
	//first search source in zgrid
	int zgrid1 = 0,zgrid2 = zgrid1+1,dz=3;  //zgird = 2
	double zbin_c[2]={0.2,0.58},mag=0;
	//double zbin_b[nzbin][2]={{0.58,0.72},{0.72,0.86},{0.86,1.02},{1.02,1.3},{1.3,2.5}};
	double zbin_b[1][2]={{0.58,2.5}};
	int ngroup;//708; 
	double (*group_pos)[4];
    group_pos = mymalloc(sizeof(double)*4*(1200*1200));
	ngroup = read_cluster(flagc,group_pos);
    printf("ngoupr=%d\n",ngroup);	

	int **found,*Nfound,pid,pindex,Nbound=0;
	double mulr=0.2; //0.2; 1.;
	double rmax,rtmp;// /4./60; //10./60
	double rmin = 50.;
	if(flagd < 1) rmax=20./60;
	else         rmax=3000;
	
	found = mymalloc(sizeof(int *)*ngroup);
	Nfound   = mymalloc(sizeof(int)*ngroup);
	memset(Nfound,0,sizeof(int)*ngroup);

	omp_set_num_threads(30);
	//for(i=0;i<1;i++)
//#pragma omp parallel for schedule(dynamic,1)
#pragma omp parallel for private(i,rtmp)
	for(i=0;i<ngroup;i++) 
	{
		if(i %10000 == 0) printf("i = %d \n",i);
		if(flagd==0) rtmp = rmax/cos(au*group_pos[i][1]);
		else rtmp = (rmax/group_pos[i][3])/au/cos(au*group_pos[i][1]);
		found[i] = linklist_search_sub_pos2(Lens,group_pos[i],rtmp,
				&Nfound[i],hoc,ll,xrange,step,NDIV,nfit);
		//printf("linklist,Nfound[%d]=%d \n",i,Nfound[i]); fflush(stdout);
	}
	printf("end \n");
	

	//------mean shear------//
	//2d:
	double rr = rmax/sqrt(2.),gstep[2];
	int ntr[2] = {20,20};
	//int ntr[2] = {10,6};
    gstep[0] = 2*rr/(ntr[0]);
    gstep[1] = gstep[0];
    SHEAR ***shear;
	shear = mymalloc(sizeof(SHEAR **)*nzbin);
	for(i=0;i<nzbin;i++)
	{
		shear[i]=shear_calloc(ntr[0],ntr[1]);
		check_shear(shear[i],ntr[0],ntr[1]);
	}
	SHEARR **shearr;
	shearr = shearr_calloc(nzbin,ntr[1]);
	
//#pragma omp parallel for private(i,j)
	for(j=0;j<nzbin;j++)
	{
		//for(i=30;i<36;i++)
		for(i=0;i<ngroup;i++)
			if(group_pos[i][2] >= zbin_c[0] && group_pos[i][2] < zbin_c[1])
			{
				if(i == 5000) printf("%d \n",i);
				mean_shear(group_pos[i],Lens,Nfound[i],found[i],shear[j],gstep,ntr,zbin_b[j],rr,rmin,flagd,flag_lg);
			}
	}
	gain_fileo_name(flagd,flagc,tmp1,tmp2);
    for(i=0;i<nzbin;i++)
	{
        //sprintf(file,"data-2017/onetime/4field-jzfield-%s-%lg-%lg-%d-%d-shear-cali2-%scra-%lg-%lg.dat",
        sprintf(file,"data-2017/onetime/w1-jzfield-%s-%lg-%lg-%d-%d-shear-cali2-%scra-%lg-%lg.dat",
                tmp2,zbin_c[0],zbin_c[1],ntr[0],ntr[1],tmp1,zbin_b[i][0],zbin_b[i][1]);
        shear_output(ntr,shear[i],shearr[i],file,gstep,rr);
    }

	return 0;
}


void gain_gstep(double rr,double rmin,double gstep[2],int ntr[2],int flag_lg)
{
    gstep[0] = 2*pi/(ntr[0]);
    if(flag_lg==0) gstep[1] = rr/ntr[1];
    else gstep[1] = (log10(rr)-log10(rmin))/ntr[1];
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

void gain_fileo_name(int flagd,int flagc,char *tmp1,char *tmp2)
{
    if(flagd==0)  sprintf(tmp1,"xy"); else sprintf(tmp1,"adxy");
    if(flagc==0) sprintf(tmp2,"clusters");
    else if(flagc==1) sprintf(tmp2,"redmapper");
    else if(flagc==-1) sprintf(tmp2,"void_mag-20_zmin0.2-rcut3");
    else if(flagc==-11) sprintf(tmp2,"void_random");
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


LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
		double xmax[2],double ymax[2])
{
	int i,n,tmp = 0;
	FILE *fp;          
	char chararr[512];
	LENS *Lens;
	myfopen(fp,file,"r");                                                               
	//fscanf(fp,"%*[^\n]%*c");  //skip the first line                                                          
	n = 0;                                                                            
	while(n < nobj)                                                                
	//while(n < 5)                                                                
	{                                                                                                                                            
		for(i=0;i<3;i++) fscanf(fp,"%s",chararr);	
		for(i=0;i<num;i++)                                                               
			fscanf(fp,"%lf",&arr[n][i]);                                                    
		//for(i=0;i<num;i++)                                                               
		//	printf("%lf ",arr[n][i]);  
		//printf("\n");   		                                          
		//if(arr[n][5] == 0 && (arr[n][2]!=0 && arr[n][3]!=0) && arr[n][6] <= 1.2) tmp++;
		if(arr[n][7] == 0 && (arr[n][4]!=0 && arr[n][5]!=0) && 
				(arr[n][12] > -99 && arr[n][12] < 99)) tmp++;
		//if(arr[n][7] == 0 && arr[n][8] <= 1.2 && arr[n][12] > -99) tmp++;
		n ++;		                                                                  
	}                                                                  
	fclose(fp);                    
	*nfit = tmp;

	printf("number of good data is %d \n",tmp);
	Lens = mymalloc(sizeof(LENS)*tmp);
	good_data(arr,nobj,tmp,Lens,zmax,xmax,ymax);
	return Lens;
}

void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2])
{
	int i,tmp = 0;
	float ymax = -10,xmax = 0,zmax = 0,xmin = 1000,ymin = 1000;
	for(i=0;i<nobj;i++)
		//if(arr[i][5] == 0 && (arr[i][2]!=0 && arr[i][3]!=0) && arr[i][6] <= 1.2)
		if(arr[i][7] == 0 && (arr[i][4]!=0 && arr[i][5]!=0) &&
				(arr[i][12] > -99 && arr[i][12] < 99))
		//if(arr[i][7] == 0 && arr[i][8] <= 1.2 && arr[i][12] > -99)
		{
			Lens[tmp].e1 = arr[i][4];  
			Lens[tmp].e2 = arr[i][5]-arr[i][10];
			Lens[tmp].m  = arr[i][9];
			Lens[tmp].mag = arr[i][11];  //absolute L
			//Lens[tmp].mag = arr[i][12];  //mag
			Lens[tmp].wt = arr[i][6]*0.01;
			Lens[tmp].pos[0] = -arr[i][0];
			Lens[tmp].pos[1] = arr[i][1];
			Lens[tmp].pos[2] = arr[i][8];
			
			if(xmax < arr[i][0]) xmax = arr[i][0];
			if(ymax < arr[i][1]) ymax = arr[i][1];
			if(zmax < arr[i][8]) zmax = arr[i][8];
			if(xmin > arr[i][0]) xmin = arr[i][0];
			if(ymin > arr[i][1]) ymin = arr[i][1];
						
			tmp ++;
		}
	printf("xmax = %f ,xmin = %f ,ymax = %f ,ymin = %f ,zmax = %f \n",xmax,xmin,ymax,ymin,zmax);
	*z = zmax;
	xma[1] = xmax+0.1;  xma[0] = xmin-0.1;
	yma[1] = ymax+0.1;  yma[0] = ymin-0.1;
}

int if_range(double xrange[2][2],double tmp[2])
{
	if(-tmp[0]>xrange[0][0] && -tmp[0]<xrange[0][1]
			&& tmp[1]>xrange[1][0] && tmp[1]<xrange[1][1])
		return 1;
	else return 0;
	
}


////select the infront galaxies
//int clusters(double (*pos)[4],int num,double xrange[2][2])
//{
//    int i,j,k,index,id,ncol,ngroup,ntmp;
//    double tmp1[2],tmp2,tmp[3];
//    FILE *fp;
//    char dir[512]=".";
//    char file[512];
//    //sprintf(file,"%s/cfht_iso_158.dat",dir);
//    sprintf(file,"%s/clusters",dir);
//    myfopen(fp,file,"r");
//
//    fscanf(fp,"%*[^\n]%*c");
//    ngroup = num;
//    for(i=0;i<num;i++)
//    {
//        fscanf(fp,"%lf %lf %lf",&tmp[0],&tmp[1],&tmp[2]);
//        pos[i][0] = -tmp[0];
//        pos[i][1] = tmp[1];
//        pos[i][2] = tmp[2];
//        pos[i][3] = DH*qromb(angu_dia_d,0,tmp[2])/(1.+tmp[2]);
//        fscanf(fp,"%lf",&tmp2);
//        fscanf(fp,"%d",&ntmp);
//        printf("%lf %lf %lf %lf %lf %d\n",pos[i][0],pos[i][1],pos[i][2],pos[i][3],tmp2,ntmp);
//    }
//    fclose(fp);
//    printf("num=%d,ngroup=%d\n",num,ngroup);
//    return ngroup;
//}


double sin2a_d(double x,double y)
{
	return 2*x*y/(x*x+y*y);
}

double cos2a_d(double x,double y)
{
	return (x*x-y*y)/(x*x+y*y);
}


void mean_shear(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
		double gstep[2],int ntr[2],double zbin[2],double rr,double rmin,int flagd,int flag_lg)
{
	int i,j,pid,grid[2],flag;
	double e1,e2,etmp[2],wt,wt2,sin2a,cos2a;
	double dx,dy,dr,theta,tmp,alpha;
#define POS(i,j) Lens[i].pos[j]
	//printf("center x,y=%f %f ,gstep[0],gstep[1]=%f %f \n",cen[0],cen[1],gstep[0],gstep[1]);
	for(i=0;i<Nfound;i++)
	{
		pid = found[i];// printf("pid=found[%d]=%d,",i,pid); fflush(stdout);
		//if(pid != cenid)
		//if(POS(pid,2) >= zbin[0] && POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1)                     //change in 20150615
		if(POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1)                     //change in 20150615
		{
			if(flagd==0)
			{
				dx = (POS(pid,0) - cen[0])*cos(au*cen[1]);
				dy = POS(pid,1) - cen[1];
			}
			else
			{
				dx = (POS(pid,0) - cen[0])*cen[3]*au*cos(au*cen[1]);
				dy = (POS(pid,1) - cen[1])*cen[3]*au;
			}
            grid[0] = floor(dx/gstep[0]) + ntr[0]/2;
            grid[1] = floor(dy/gstep[1]) + ntr[1]/2;
            if(fabs(dx) < rr && fabs(dy) < rr)
            {
                e1 = Lens[pid].e1;
                e2 = Lens[pid].e2;
                wt = Lens[pid].wt;
                wt2 = wt*wt;
                shear[grid[0]][grid[1]].gt +=e1*wt;  //here et=e1
                shear[grid[0]][grid[1]].gr +=e2*wt;  //here er=e2   
                shear[grid[0]][grid[1]].wt +=wt;
                shear[grid[0]][grid[1]].m  +=(1+Lens[pid].m)*wt;
                shear[grid[0]][grid[1]].nrr ++;
                shear[grid[0]][grid[1]].err1 +=e1*e1*wt2;
                shear[grid[0]][grid[1]].err2 +=e2*e2*wt2;
                shear[grid[0]][grid[1]].dr   += dx*wt;         //here dr=dx
                shear[grid[0]][grid[1]].angle += dy*wt;        //here angle=dy
            }               
        }
    }

#undef POS
}


double et1(double e1,double e2,double sin2a,double cos2a)
{
	return -e1*cos2a-e2*sin2a;
}

double et2(double e1,double e2,double sin2a,double cos2a)
{
	return e1*sin2a-e2*cos2a;
}


void check1d(int *arr,int n)
{
	int i;
	for(i=0;i<n;i++)
		if(arr[i] != 0) printf("use of memset is wrong!\n");
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

double gain_r_grid(int flag_lg,double gstep1,double rmin,int j)
{
    if(flag_lg == 0) return gstep1*(j+0.5);
    else return pow(10.,log10(rmin)+gstep1*(j+0.5));
}

void shear_output(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],double rr)
{
    int i,j;
    FILE *fp;

    myfopen(fp,file,"w");

    for(i=0;i<ntr[0];i++)
        for(j=0;j<ntr[1];j++)
        {
            if(shear[i][j].nrr > 0)
            {
                shear[i][j].gt /=shear[i][j].m;
                shear[i][j].gr /=shear[i][j].m;
                shear[i][j].dr /=shear[i][j].wt;
                shear[i][j].angle /=shear[i][j].wt;
                shear[i][j].err1 = sqrt(shear[i][j].err1);
                shear[i][j].err2 = sqrt(shear[i][j].err2);
                shear[i][j].err1 /=shear[i][j].wt;
                shear[i][j].err2 /=shear[i][j].wt;
            }
            fprintf(fp,"%f %f %f %f %d %f %f %f %f\n",
                    shear[i][j].gt,shear[i][j].gr,gstep[0]*(i+0.5)-rr,gstep[1]*(j+0.5)-rr,
                    shear[i][j].nrr,shear[i][j].err1,shear[i][j].err2,shear[i][j].dr,
                    shear[i][j].angle);
        }

    fclose(fp);

}

#undef omegal0 
#undef omegam0 
#undef H0 
#undef DH 
#undef au 
