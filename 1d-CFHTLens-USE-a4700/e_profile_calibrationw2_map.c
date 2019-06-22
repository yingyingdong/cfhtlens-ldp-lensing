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

double **gain_xgrid(int NDIV[3],double xrange[2][2]);
void write_grid_num(SHEAR ***shear,int ngroup);
int *gain_grid_flag(char *file,int nx,int ny);
void gain_xrange_here(int field,double xrange[2][2]);
void gain_group_pos(int NDIV[3],double **xgrid,double (*group_pos)[4]);
void get_line(char *file,int *nobj);
LENS *load_data(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
		        double xmax[2],double ymax[2]);
double **D_malloc2_here(int n,int NDIV[2]);
void good_data(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2]);
//int clusters(double (*pos)[4],int num,double xrange[2][2]);
float arcsin(double x,double y);
void mean_shear(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
		double gstep[2],int ntr[2],double zbin[2],double rmin,int flagd);
double et1(double e1,double e2,double sin2a,double cos2a);
double et2(double e1,double e2,double sin2a,double cos2a);
static double step[2];
void check1d(int *arr,int n);
void check_shear(SHEAR **shear,int n1,int n2);
void shear_output(int ntr[2],SHEAR ***shear,SHEARR **shearr,char *file,double gstep[2],
        double rmin,int ngroup,double (*group_pos)[4]);

#define omegal0 0.7
#define omegam0 0.3
#define H0 72.
#define au (pi/180.)  //error
#define nzbin 5  //error
//double angu_dia_d(double z)
//{
//	return 1./sqrt(omegal0+omegam0*(1+z)*(1+z)*(1+z));
//}


int main(int argc,char **argv)
{

	FILE *fp;
	//char file[1024]="../../22/hebing/all.tsv";
	char dir[1024]="/home/fydong/data_dir/data/CFHTLens/";
    char field[4][1024]={"CFHTLens_2015-07-15T01:50:34.tsv","CFHTLens_2015-07-22T03_52_26.tsv",
        "CFHTLens_2015-07-14T02:12:29.tsv","CFHTLens_2015-12-18T03_22_58.tsv"};
	char tmp1[512],tmp2[512],tmp3[512],buf[1024],name[512],file[1024];
    char file_flag[512]="/home/fydong/data_dir/data/mask_grid/100.dat",
         ffname[4][16]={"w1","w2","w3","w4"};
	int i,j,nfit,flagc=1,flagd=0,nfield=2;
    double zbin[2]={0.1,2.};
	
	
	//------load data------//
	int nobj = 0,num=13;
	int NDIV[2]={200,200},*grid_flag;
    double **arr;
	double zmax,xmax[2],ymax[2];
	LENS *Lens;
    sprintf(file,"%s%s",dir,field[nfield]);
	get_line(file,&nobj);
	//nobj = 22724706;
	arr = D_malloc2(nobj,num);	
	Lens = load_data(file,nobj,num,arr,&nfit,&zmax,xmax,ymax);	
	free_D_malloc2(arr,nobj);
    //grid_flag = gain_grid_flag(file_flag,NDIV[0],NDIV[1]);
	printf("zmax = %f \n",zmax);



	//------make grid------//
	double xrange[2][2],**xgrid;
    int ngroup,ngroup_tmp; double (*group_pos)[4];

	xgrid = D_malloc2_here(2,NDIV);
	for(i=0;i<2;i++)
		xrange[1][i] = ymax[i];
	xrange[0][0] = -xmax[1];  
	xrange[0][1] = -xmax[0];
    gain_xrange_here(nfield,xrange);
    xgrid = gain_xgrid(NDIV,xrange);

	printf("xmax:%f %f ,ymax:%f %f ,zmax:%f\n",xrange[0][0],xrange[0][1],
			xrange[1][0],xrange[1][1],zmax);
    ngroup = NDIV[0]*NDIV[1];
    group_pos = mymalloc(sizeof(double)*4*ngroup);
    gain_group_pos(NDIV,xgrid,group_pos);
	
	//------make linklist------//	
	int *hoc,*ll;
	hoc = mymalloc(sizeof(int)*NDIV[0]*NDIV[1]);
	ll  = mymalloc(sizeof(int)*nfit);
	
	for(i=0;i<2;i++) step[i] = (xrange[i][1]-xrange[i][0])/NDIV[i];
	makell_sub(Lens,nfit,hoc,ll,NDIV,xrange,step,xgrid);

	//check_num(Lens,NDIV,xgrid,hoc,ll);
		
	//search grid
	//first search source in zgrid
	int **found,*Nfound,pid,pindex,Nbound=0;
	double rmax,rtmp;// /4./60; //10./60
	double rmin = 0.;//0.1/60;
    rmax = 20./60;//MAX_D(step[0],step[1]);//5./60;

	found = mymalloc(sizeof(int *)*ngroup);
	Nfound   = mymalloc(sizeof(int)*ngroup);
	memset(Nfound,0,sizeof(int)*ngroup);

    double rr = rmax,gstep[2];
    int ntr[2] = {1,1};
    gstep[0] = 2*pi/(ntr[0]);
    gstep[1] = rr/ntr[1];
    SHEAR ***shear; SHEARR **shearr;
    shear  = mymalloc(sizeof(SHEAR **)*ngroup);
    shearr = shearr_calloc(ngroup,ntr[1]); 

    for(i=0;i<ngroup;i++)
    {
        shear[i]=shear_calloc(ntr[0],ntr[1]);
        check_shear(shear[i],ntr[0],ntr[1]);
    }

	omp_set_num_threads(30);
//#pragma omp parallel for schedule(dynamic,1)
#pragma omp parallel for private(i,rtmp)
	for(i=0;i<ngroup;i++)
    //if(grid_flag[i]==1)
	{
		if(i %10000 == 0) printf("i = %d \n",i);
		if(flagd==0) rtmp = rmax/cos(au*group_pos[i][1]);
		else rtmp = (rmax/group_pos[i][3])/au/cos(au*group_pos[i][1]);
		found[i] = linklist_search_sub_pos2(Lens,group_pos[i],rtmp,
				&Nfound[i],hoc,ll,xrange,step,NDIV,nfit);
        if(Nfound[i]>0)
            mean_shear(group_pos[i],Lens,Nfound[i],found[i],shear[i],gstep,ntr,zbin,rmin,flagd);
	}
	printf("end \n");
	

	

	if(flagd==0)  sprintf(tmp1,"r");
	else          sprintf(tmp1,"adr");
    //sprintf(file,"data-2017/EB_map/w1-zbin%lg-%lg-%d-%d-%s-rmax%lg-rmin%lg-shear-cali2-cra-nomaglimit.dat",
    sprintf(file,"data-2017/EB_map/%s-%d-%d-zbin%lg-%lg-%d-%d-%s-rmax%lg-rmin%lg-shear-cali2-cra-nogridflag-nomaglimit.dat",
            ffname[nfield],NDIV[0],NDIV[1],zbin[0],zbin[1],ntr[0],ntr[1],tmp1,rmax,rmin); // 05 means : gird 0-5 
    shear_output(ntr,shear,shearr,file,gstep,rmin,ngroup,group_pos);	

	return 0;
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

int *gain_grid_flag(char *file,int nx,int ny)
{
    FILE *fp;
    int i,j,*flag,n=nx*ny;
    flag = mymalloc(sizeof(int)*n);

    myfopen(fp,file,"r");
    for(i=0;i<n;i++)
        fscanf(fp,"%d",&flag[i]);
    fclose(fp);
    return flag;
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
    myfopen(fp,"group_pos.dat","w");
    for(i=0;i<NDIV[1]*NDIV[0];i++)
    {
        fprintf(fp,"%lg %lg\n",group_pos[i][0],group_pos[i][1]);
    }
    fclose(fp);
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
		double xmax[2],double ymax[2])
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
		//if(arr[n][5] == 0 && (arr[n][2]!=0 && arr[n][3]!=0) && arr[n][6] <= 1.2) tmp++;
		if(arr[n][7] == 0 && (arr[n][4]!=0 && arr[n][5]!=0)) tmp++;
		//	&&	(arr[n][12] > -99 && arr[n][12] < 99)) tmp++;
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
		if(arr[i][7] == 0 && (arr[i][4]!=0 && arr[i][5]!=0))
		//		&& (arr[i][12] > -99 && arr[i][12] < 99))
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
		double gstep[2],int ntr[2],double zbin[2],double rmin,int flagd)
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
		if(POS(pid,2) >= zbin[0] && POS(pid,2)<zbin[1])                     //change in 20150615
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
			dr = sqrt(dx*dx + dy*dy);
			e1 = Lens[pid].e1;
			e2 = Lens[pid].e2;
			theta = atan2(dy,dx);
			if(theta < 0) tmp = 2*pi+theta;
			else tmp = theta;
			grid[0] = floor(tmp/gstep[0]); 
			//grid[1] = floor( (log10(dr)-log10(rmin))/gstep[1] );  
			grid[1] = floor(dr/gstep[1]);
			sin2a = sin2a_d(dx,dy); 
			cos2a = cos2a_d(dx,dy);
            wt = Lens[pid].wt;
			if(dr>=rmin && grid[1] >= 0 && grid[1]<ntr[1] && Lens[pid].wt>0)  
			{
				Lens[pid].angle = tmp;
				etmp[0] = et1(e1,e2,sin2a,cos2a);
				etmp[1] = et2(e1,e2,sin2a,cos2a);
				wt2 = wt*wt;
				shear[grid[0]][grid[1]].gt +=etmp[0]*wt;
				shear[grid[0]][grid[1]].gr +=etmp[1]*wt;
				shear[grid[0]][grid[1]].wt +=wt;
				shear[grid[0]][grid[1]].m  +=(1+Lens[pid].m)*wt;
				shear[grid[0]][grid[1]].nrr ++;
				shear[grid[0]][grid[1]].err1 +=etmp[0]*etmp[0]*wt2;
				shear[grid[0]][grid[1]].err2 +=etmp[1]*etmp[1]*wt2;
				shear[grid[0]][grid[1]].dr   +=dr*wt;
				shear[grid[0]][grid[1]].angle +=tmp*wt;
				//printf("dx=%f,dy=%f,x=%f,y=%f,theta=%f,grid[0]=%d,grid[1]=%d,gridi=%f,gridj=%f,e1,e2=%f %f,et,er=%f %f \n",
				//		dx,dy,POS(pid,0),POS(pid,1),tmp,grid[0],grid[1],gstep[0]*grid[0],gstep[1]*grid[1],e1,e2,et1(e1,e2,alpha),et2(e1,e2,alpha));
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

void shear_output(int ntr[2],SHEAR ***shear,SHEARR **shearr,char *file,double gstep[2],double rmin,
        int ngroup,double (*group_pos)[4])
{   
	int i,j,k;
	FILE *fp;

	myfopen(fp,file,"w");
/*	for(i=0;i<ntr[1];i++)
		for(j=0;j<ntr[0];j++)
		{
			shearr[i].gt += shear[j][i].gt;
			shearr[i].gr += shear[j][i].gr;
			shearr[i].dr += shear[j][i].dr;
			shearr[i].m  += shear[j][i].m;
			shearr[i].err1 += shear[j][i].err1;
			shearr[i].err2 += shear[j][i].err2;
			shearr[i].wt += shear[j][i].wt;
			shearr[i].nrr += shear[j][i].nrr;
		}
	for(i=0;i<ntr[1];i++) if(shearr[i].wt>0)
	{
		shearr[i].gt /=shearr[i].m;
		shearr[i].gr /=shearr[i].m;
		shearr[i].dr /=shearr[i].wt;
		shearr[i].err1 =sqrt(shearr[i].err1);
		shearr[i].err2 =sqrt(shearr[i].err2);
		shearr[i].err1 /=shearr[i].wt;
		shearr[i].err2 /=shearr[i].wt;
	}

*/
    for(k=0;k<ngroup;k++)
        for(i=0;i<ntr[0];i++)
	    	for(j=0;j<ntr[1];j++)
	    	{
	    		if(shear[k][i][j].nrr > 0)
	    		{
	    			shear[k][i][j].gt /=shear[k][i][j].m;
	    			shear[k][i][j].gr /=shear[k][i][j].m;
	    			shear[k][i][j].dr /=shear[k][i][j].wt;
	    			shear[k][i][j].angle /=shear[k][i][j].wt;
	    			shear[k][i][j].err1 = sqrt(shear[k][i][j].err1);
	    			shear[k][i][j].err2 = sqrt(shear[k][i][j].err2);
	    			shear[k][i][j].err1 /=shear[k][i][j].wt;
	    			shear[k][i][j].err2 /=shear[k][i][j].wt;
	    		}
	    		fprintf(fp,"%f %f %f %f %d %f %f %f %f\n",
	    		shear[k][i][j].gt,shear[k][i][j].gr,group_pos[k][0],group_pos[k][1],
	    		shear[k][i][j].nrr,shear[k][i][j].err1,shear[k][i][j].err2,shear[k][i][j].angle,
	    		shear[k][i][j].dr);
	        }

	fclose(fp);
/*	
	for(i=0;i<ntr[1];i++)
	printf("%f %f %f %d %f %f %f\n",
		shearr[i].gt,shearr[i].gr,shearr[i].wt,shearr[i].nrr,shearr[i].dr,shearr[i].err1,shearr[i].err2);
*/
}

#undef omegal0 
#undef omegam0 
#undef H0 
#undef DH 
#undef au 
