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
#include "fft_m.h"
#include "get_z_dal.h"
#include "thread_assign.h"


void good_data_o(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2])
{
	int i,tmp = 0;
	float ymax = -10,xmax = 0,zmax = 0,xmin = 1000,ymin = 1000;
	for(i=0;i<nobj;i++)
		//if(arr[i][5] == 0 && (arr[i][2]!=0 && arr[i][3]!=0) && arr[i][6] <= 1.2)
		if(arr[i][7] == 0 && (arr[i][4]!=0 && arr[i][5]!=0) &&
				(arr[i][12] > -99 && arr[i][12] < 99))
		//if(arr[i][7] == 0 && arr[i][8] <= 1.2 && arr[i][12] > -99)
		{
			Lens[tmp].pz = arr[i][14]; 
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


LENS *load_data_o(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
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
	//while(n < 100)                                                                
	{                                                                                                                                            
		for(i=0;i<3;i++) fscanf(fp,"%s",chararr);	
		for(i=0;i<num;i++)                                                               
			fscanf(fp,"%lg",&arr[n][i]);                                                    
		//for(i=0;i<num;i++)                                                               
		//	printf("%lg ",arr[n][i]); printf("\n"); 
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
	good_data_o(arr,nobj,tmp,Lens,zmax,xmax,ymax);
	return Lens;
}


double et1o(double e1,double e2,double sin2a,double cos2a)
{
	return -e1*cos2a-e2*sin2a;
}

double et2o(double e1,double e2,double sin2a,double cos2a)
{
	return e1*sin2a-e2*cos2a;
}


double sin2a_do(double x,double y)
{
	return 2*x*y/(x*x+y*y);
}

double cos2a_do(double x,double y)
{
	return (x*x-y*y)/(x*x+y*y);
}


void mean_shear_r(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
		double gstep[2],int ntr[2],double zbin[2],double rmin, 
       double *DAzlg,double zmean,double zsteplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density)
{
	int i,j,pid,grid[2],flag;
	double e1,e2,etmp[2],wt,wt2,sin2a,cos2a;
	double dx,dy,dr,theta,tmp,alpha,sigmaz,a=1,density;

    if(flag_com == 1) a=1./(1+cen[2]);
#define POS(i,j) Lens[i].pos[j]
	//printf("center x,y=%f %f ,gstep[0],gstep[1]=%f %f \n",cen[0],cen[1],gstep[0],gstep[1]);
	for(i=0;i<Nfound;i++)
	{
		pid = found[i];// printf("pid=found[%d]=%d,",i,pid); fflush(stdout);
		//if(pid != cenid)
		//if(POS(pid,2) >= zbin[0] && POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1 && Lens[pid].pz>0 && Lens[pid].pz<pz)                     //change in 20150615
		if(POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1 && Lens[pid].pz>0 && Lens[pid].pz<pz)                     //change in 20150615
		//if(POS(pid,2)<zbin[1] && POS(pid,2)>cen[2] && Lens[pid].pz>0 && Lens[pid].pz<pz)                     //change in 20150615
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
			dr = sqrt(dx*dx + dy*dy)/a;
			e1 = Lens[pid].e1;
			e2 = Lens[pid].e2;
			theta = atan2(dy,dx);
			if(theta < 0) tmp = 2*pi+theta;
			else tmp = theta;
			grid[0] = floor(tmp/gstep[0]); 
			if(flag_lg == 0) grid[1] = floor(dr/gstep[1]);
            else grid[1] = floor( (log10(dr)-log10(rmin))/gstep[1] );
			sin2a = sin2a_do(dx,dy); 
			cos2a = cos2a_do(dx,dy);
            wt = Lens[pid].wt;
			if(wt>0 && dr>=rmin && grid[1] >= 0 && grid[1]<ntr[1])  
			{
				sigmaz = sigma(cen[2],POS(pid,2),zmean,zsteplg,cen[3],DAzlg)*a*a;
                if(flag_wt==1) wt = wt/(sigmaz*sigmaz);
                if(flag_density == 1) density=sigmaz/1e12; else density=1.;
                Lens[pid].angle = tmp;
				etmp[0] = et1o(e1,e2,sin2a,cos2a)*wt*density;
				etmp[1] = et2o(e1,e2,sin2a,cos2a)*wt*density;
				wt2 = wt*wt;
                shear[grid[0]][grid[1]].gt +=etmp[0];
				shear[grid[0]][grid[1]].gr +=etmp[1];
				shear[grid[0]][grid[1]].wt +=wt;
				shear[grid[0]][grid[1]].m  +=(1+Lens[pid].m)*wt;
				shear[grid[0]][grid[1]].nrr ++;
				shear[grid[0]][grid[1]].err1 +=etmp[0]*etmp[0];
				shear[grid[0]][grid[1]].err2 +=etmp[1]*etmp[1];
				shear[grid[0]][grid[1]].dr   +=dr*wt;
				shear[grid[0]][grid[1]].angle +=tmp*wt;
			}
		}
	}

#undef POS
}

void mean_shear_xy(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
		double gstep[2],int ntr[2],double zbin[2],double rr,double rmin,
        double *DAzlg,double zmin,double zsteplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density)
{
	int i,j,pid,grid[2],flag;
	double e1,e2,etmp[2],wt,wt2,sin2a,cos2a;
	double dx,dy,dr,theta,tmp,alpha,a=1,density,sigmaz;

    if(flag_com == 1) a=1./(1+cen[2]);
#define POS(i,j) Lens[i].pos[j]
	//printf("center x,y=%f %f ,gstep[0],gstep[1]=%f %f \n",cen[0],cen[1],gstep[0],gstep[1]);
	for(i=0;i<Nfound;i++)
	{
		pid = found[i];// printf("pid=found[%d]=%d,",i,pid); fflush(stdout);
		//if(pid != cenid)
		//if(POS(pid,2) >= zbin[0] && POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1)                     //change in 20150615
		if(POS(pid,2)<zbin[1] && POS(pid,2)>cen[2]+0.1  && Lens[pid].pz>0 && Lens[pid].pz<pz)                     //change in 20150615
		{
			if(flagd==0)
			{
				dx = (POS(pid,0) - cen[0])*cos(au*cen[1])/a;
				dy = (POS(pid,1) - cen[1])/a;
			}
			else
			{
				dx = (POS(pid,0) - cen[0])*cen[3]*au*cos(au*cen[1])/a;
				dy = (POS(pid,1) - cen[1])*cen[3]*au/a;
			}
            grid[0] = floor(dx/gstep[0]) + ntr[0]/2;
            grid[1] = floor(dy/gstep[1]) + ntr[1]/2;
            if(fabs(dx) < rr && fabs(dy) < rr)
            {
                wt = Lens[pid].wt;
                wt2 = wt*wt;
                sigmaz = sigma(cen[2],POS(pid,2),zmin,zsteplg,cen[3],DAzlg)*a*a;
                if(flag_wt==1) wt = wt/(sigmaz*sigmaz);
                if(flag_density == 1) density=sigmaz/1e12; else density=1.;
                e1 = Lens[pid].e1*density*wt;
                e2 = Lens[pid].e2*density*wt;
                shear[grid[0]][grid[1]].gt +=e1;  //here et=e1
                shear[grid[0]][grid[1]].gr +=e2;  //here er=e2   
                shear[grid[0]][grid[1]].wt +=wt;
                shear[grid[0]][grid[1]].m  +=(1+Lens[pid].m)*wt;
                shear[grid[0]][grid[1]].nrr ++;
                shear[grid[0]][grid[1]].err1 +=e1*e1;
                shear[grid[0]][grid[1]].err2 +=e2*e2;
                shear[grid[0]][grid[1]].dr   += dx*wt;         //here dr=dx
                shear[grid[0]][grid[1]].angle += dy*wt;        //here angle=dy
            }               
        }
    }

#undef POS
}


void gain_gstep_o(double rr,double rmin,double gstep[2],int ntr[2],int flag_lg)
{
    gstep[0] = 2*pi/(ntr[0]);
    if(flag_lg==0) gstep[1] = rr/ntr[1];
    else gstep[1] = (log10(rr)-log10(rmin))/ntr[1];
}

void gain_nr_gstep_o(int ntr[2],double gstep[2],int flag_lg,int flagd,int flag_xy,double *rmin,double *rmax)
{
    if(flag_xy == 0)
    {
        ntr[0]=1; ntr[1]=60;
        if(flagd < 1) ntr[1] = 20; if(flag_lg==1) ntr[1]=20;
        if(flagd < 1) {*rmax=20./60;*rmin=0.1/10;}
        else          {*rmax=12000; *rmin=10;}
        gain_gstep_o(*rmax,*rmin,gstep,ntr,flag_lg);
    }
    else
    {
        ntr[0]=ntr[1]=20;
        if(flagd < 1) *rmax=70./60;
        else *rmax = 12000.;
        gstep[0] = 2*(*rmax)/sqrt(2.)/(ntr[0]);
        gstep[1] = gstep[0];
    }
}

double gain_r_grid_o(int flag_lg,double gstep1,double rmin,int j)
{
    if(flag_lg == 0) return gstep1*(j+0.5);
    else return pow(10.,log10(rmin)+gstep1*(j+0.5));
}

void shear_output_r(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],
        double rmin,int flag_lg,int nuse)
{   
	int i,j; double tmp;
	FILE *fp;

    myfopen(fp,file,"w");
    for(i=0;i<ntr[1];i++)
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
			tmp = gain_r_grid_o(flag_lg,gstep[1],rmin,j);
            fprintf(fp,"%lg %lg %lg %lg %d %lg %lg %lg %lg %e\n",
			shear[i][j].gt,shear[i][j].gr,gstep[0]*(i+0.5),tmp,
			shear[i][j].nrr,shear[i][j].err1,shear[i][j].err2,shear[i][j].angle,
			shear[i][j].dr,shear[i][j].wt/nuse);
	    }

	fclose(fp);
	
	for(i=0;i<ntr[1];i++)
	printf("%lg %lg %lg %d %lg %lg %lg\n",
		shearr[i].gt,shearr[i].gr,shearr[i].wt,shearr[i].nrr,shearr[i].dr,shearr[i].err1,shearr[i].err2);
}

void shear_output_xy(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],double rr)
{   
	int i,j; double tmp;
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
            fprintf(fp,"%lg %lg %lg %lg %d %lg %lg %lg %lg\n",
                    shear[i][j].gt,shear[i][j].gr,gstep[0]*(i+0.5)-rr,gstep[1]*(j+0.5)-rr,
                    shear[i][j].nrr,shear[i][j].err1,shear[i][j].err2,shear[i][j].dr,
                    shear[i][j].angle);
        }

    fclose(fp);


}


void save_found_zbin_o(int **Nfound,int nzbin,int nx,int ny,char *name)
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

void save_void_pos_o(int **Nfound,int **flag_mask,int NDIV[2],int nzbin,double (*group_pos)[4],
        char *name,double zstep,double zmean,double xc[2])
{
    int i,j,k,nelem,num=0;
    double ra;
    FILE *fp;

    myfopen(fp,name,"w")
    for(k=0;k<nzbin;k++)
        for(i=0;i<NDIV[1];i++) //hang
        {
            for(j=0;j<NDIV[0];j++) //lie
            {
                nelem = ELEM(i,j,NDIV[0]);
                if(flag_mask[i][j] == 1 && Nfound[nelem][k] == 0)
                {
                    ra = (group_pos[nelem][0]-xc[0])/cos(group_pos[nelem][1]*au)+xc[0];
                    fprintf(fp,"%e %e %e %e\n",ra,group_pos[nelem][0],group_pos[nelem][1],zmean);
                    num++;//printf("there is!\n");
                }
            }
        }
    fclose(fp);
    printf("num=%d\n",num);
}

void save_void_pos_largegrid_o(int **Nfound,int NDIV[2],int nzbin,double (*group_pos)[4],
        char *name,double zstep,double zmean,double xc[2])
{
    int i,j,k,nelem,num=0;
    double ra;
    FILE *fp;

    myfopen(fp,name,"w")
    for(k=0;k<nzbin;k++)
        for(i=0;i<NDIV[1];i++) //hang
        {
            for(j=0;j<NDIV[0];j++) //lie
            {
                nelem = ELEM(i,j,NDIV[0]);
                //if(flag_mask[i][j] == 1 && Nfound[nelem][k] == 0)
                if(Nfound[nelem][k] == 0)
                {
                    ra = (group_pos[nelem][0]-xc[0])/cos(group_pos[nelem][1]*au)+xc[0];
                    fprintf(fp,"%e %e %e %e\n",ra,group_pos[nelem][0],group_pos[nelem][1],zmean);
                    num++;//printf("there is!\n");
                }
            }
        }
    fclose(fp);
    printf("num=%d\n",num);
}



void gain_file_name_xy(int flag_xy,int flag_lg,int flagd,int flagc,int wt_sigma,
        int flag_com,char *tmp1,char *tmp2,char *tmpwt,int ngrp_use,
        double mag_limit,double zmean_v,double rcut,int nzbin_v,double zstep)
{
    char tmp[512];

    if(flag_xy == 0)
    {
        if(flag_lg==0)
            {if(flagd==0)  sprintf(tmp,"r"); else sprintf(tmp,"adr");}
        else {if(flagd==0) sprintf(tmp,"lgr"); else sprintf(tmp,"lgadr");};
    }
    else
    {
        if(flag_lg==0)
            {if(flagd==0)  sprintf(tmp,"xy"); else sprintf(tmp,"adxy");}
        else {if(flagd==0) sprintf(tmp,"lgxy"); else sprintf(tmp,"lgadxy");};
    }

    if(flag_com == 1) sprintf(tmp1,"%s-COM-mul%d",tmp,ngrp_use);
    else sprintf(tmp,"%s-mul%d",tmp,ngrp_use);

    if(flagc==-1) sprintf(tmp2,"void_mag%lg_zmean%lg-rcut%lg-nbin%d-zstep%lg-step0.006139",mag_limit,zmean_v,rcut,nzbin_v,zstep);
    else if(flagc==-10) sprintf(tmp2,"void_mag-20_zmin0.2-rcut3-fofr1POS-n4"); //n20:at least 20 point of one fof
    //else if(flagc==-10) sprintf(tmp2,"void_mag-20_zmin0.2-rcut3-fofr0.4pos");
    else if(flagc==-11) sprintf(tmp2,"void_randommul20");
    if(wt_sigma == 0) sprintf(tmpwt,"wt"); 
    else sprintf(tmpwt,"wtsigma"); 


}

int my_void_here_o(double (*pos)[5],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg)
{
    int nline,i,nz2; double tmp;
    FILE *fp;
    //char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/w1",file[1024];
    char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/total",file[1024];

    //sprintf(file,"%s/w1_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmin%lg-rmax%lg_magl%lg-pos.dat",dir,nzbin_v,zstep,zmin_v,rcut/60.,mag_limit);
    sprintf(file,"%s/4_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmean%lg-rmax%lg_magl%lg-maskratio%lg-pos.dat",
    //sprintf(file,"%s/4_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmin%lg-rmax%lg_magl%lg-maskratio%lg-pos-randompos_mul1.dat",
            dir,nzbin_v,zstep,zmean_v,rcut/60.,mag_limit,mask_ratio);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg %lg",&pos[i][0],&tmp,&pos[i][1],&pos[i][2]); //for void
        //fscanf(fp,"%lg %lg %lg",&pos[i][0],&pos[i][1],&pos[i][2]);  //random for void
        //pos[i][2] = zmean_v;
        nz2 = floor( (log10(1+pos[i][2])-log10(1+zmin))/zsteplg );
        pos[i][3] = DAzlg[nz2];
        if(i<100)
        {   printf("%d,%lg,%lg\n",i,pos[i][3],DH*qromb(angu_dia_d,0,pos[i][2])/(1.+pos[i][2]));}
        //pos[i][3] = DH*qromb(angu_dia_d,0,pos[i][2])/(1.+pos[i][2]);
    }
    fclose(fp);
    return nline;
}


int my_void_meanfind_o(double (*pos)[5],double rcut,double zstep,double zmean_v,int nzbin_v,double mag_limit,
        double mask_ratio,double zmin,double zsteplg,double *DAzlg,double xrange[2][2],double cutratio)
{
    int nline,i,j,nz2,nfind,nn=0,nmean=0; double tmp=0,mean=0;
    FILE *fp;
    //char dir[1024] = "../1d-CFHTLens-USE-a4700/data-2017/mask_flag/maskmethod1_USE/w1",file[1024];
    char dir[1024] = "../1d-CFHTLens-USE-a4700/data-2017/mask_flag/maskmethod1_USE/total",file[1024];
    double (*g_pos)[5];
    g_pos = mymalloc(sizeof(double)*5*(2000*2000));


    sprintf(file,"%s/4_maskmethod1_nzbin%d-step0.006139-zstep%lg-zmean%lg-rmax%lg_magl%lg-maskratio%lg-nfoundpos.dat",
            dir,nzbin_v,zstep,zmean_v,rcut/60.,mag_limit,mask_ratio);
    get_line(file,&nline);
    myfopen(fp,file,"r");

    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg %lg %d",&g_pos[nn][0],&tmp,&g_pos[nn][1],&g_pos[nn][2],&nfind); //void
        g_pos[nn][4] = nfind;
        if(nfind >= 0) nn++; if(nfind>0){nmean++; mean +=nfind;}
    }
    fclose(fp);
    mean /=nmean;
    
    nmean = 0;
    for(i=0;i<nn;i++)
        if(g_pos[i][4] <= mean*cutratio)
        {
            for(j=0;j<4;j++) pos[nmean][j] = g_pos[i][j];
            nz2 = floor( (log10(1+pos[nmean][2])-log10(1+zmin))/zsteplg );
            pos[nmean][3] = DAzlg[nz2];
            nmean++;
        }

        //if(pos[nn][0]<=xrange[0][0] || pos[nn][0]>=xrange[0][1] || pos[nn][1]<=xrange[1][0] || pos[nn][1]>=xrange[1][1])
        //   printf("wrong! pos[%d][0]=%e,pos[nn][1]=%e\n",nn,pos[nn][0],pos[nn][1]); fflush(stdout);
    free(g_pos);
    return nmean;
}



int l_check(double x,double ex1,double ex2,double limit)
{
    if(x-ex1 > limit && ex2-x > limit) return 1;
    else return 0;
}

int edge_check(double x,double ex1,double ex2,double y,double ey1,double ey2,double limit)
{
    return l_check(x,ex1,ex2,limit)*l_check(y,ey1,ey2,limit);
}

            



