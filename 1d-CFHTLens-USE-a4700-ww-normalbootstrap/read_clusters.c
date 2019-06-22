#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fftw3.h"
#include "func.h"
#include "struct.h"
#include "Romberg.h"
#include "power.h"

double angu_dia_d_ww(double z)
{   
	return 1./sqrt(omegal0*pow(1.+z,3*(1.+ww))+omegam0*(1+z)*(1+z)*(1+z));
}

int cmass(double (*pos)[5],double zbin_c[2],char *fname)
{
    int i,j,k,num,ngroup;
    double tmp1[2],tmp2,tmp[3];
    FILE *fp;
    char dir[512]="/home/dfy/data/astro/clusters-caa/cmass";
    char file[512];
    sprintf(file,"%s/%s.dat",dir,fname);
    GET_line(file,&num);
    num--;
    myfopen(fp,file,"r");
    ngroup = 0;
    for(i=0;i<num;i++)
    {
        fscanf(fp,"%lf %lf %lf",&tmp[0],&tmp[1],&tmp[2]);
        pos[ngroup][0] = -tmp[0];
        pos[ngroup][1] = tmp[1];
        pos[ngroup][2] = tmp[2];
        pos[ngroup][3] = DH*qromb(angu_dia_d_ww,0,tmp[2])/(1.+tmp[2]);
        printf("%lf %lf %lf %lf\n",pos[ngroup][0],pos[ngroup][1],pos[ngroup][2],pos[ngroup][3]);
        if(pos[ngroup][2]>zbin_c[0] && pos[ngroup][2]<zbin_c[1]) ngroup++;
    }
    fclose(fp);
    printf("num=%d,ngroup=%d\n",num,ngroup);
    return ngroup;
}

int clusters(double (*pos)[5],double zbin_c[2])
{
	int i,j,k,index,id,ncol,ngroup,ntmp,num;
	double tmp1[2],tmp2,tmp[3];
	FILE *fp;
	char dir[512]="/home/dfy/data/astro/clusters-caa/CFHTLens";
	char file[512];
	sprintf(file,"%s/clusters",dir);
    GET_line(file,&num);
    num--;
    myfopen(fp,file,"r");
	fscanf(fp,"%*[^\n]%*c");
	ngroup = 0;
    for(i=0;i<num;i++)
	{
		fscanf(fp,"%lf %lf %lf",&tmp[0],&tmp[1],&tmp[2]);
		pos[ngroup][0] = -tmp[0];
		pos[ngroup][1] = tmp[1];
		pos[ngroup][2] = tmp[2];
		pos[ngroup][3] = DH*qromb(angu_dia_d_ww,0,tmp[2])/(1.+tmp[2]);
		fscanf(fp,"%lf",&tmp2);
        pos[ngroup][4] = log10(pow(10.,(0.161)*tmp2+12.39));
		fscanf(fp,"%d",&ntmp);
		printf("%lf %lf %lf %lf %lf %lf %d\n",pos[ngroup][0],pos[ngroup][1],pos[ngroup][2],pos[ngroup][3],pos[ngroup][4],tmp2,ntmp);
	    if(pos[ngroup][2]>zbin_c[0] && pos[ngroup][2]<zbin_c[1]) ngroup++;
    }
	fclose(fp);
	printf("num=%d,ngroup=%d\n",num,ngroup);
	return ngroup;
}
int redmapper(double (*pos)[5],double zbin_c[2])
{
	int i,j,k,index,id,ncol,ngroup,ntmp,num;
	double tmp1[2],tmp2,tmp[3],ntmp1;
	FILE *fp;
	char dir[512]="/home/dfy/data/astro/clusters-caa/redmapper";
	char file[512];

	//sprintf(file,"%s/cfht_iso_158.dat",dir);
	//sprintf(file,"%s/redmapper_dr8_public_v5.10_catalog.dat",dir);
	sprintf(file,"%s/4wredmapper_dr8_public_v5.10_catalog.dat",dir);
	GET_line(file,&num);
    myfopen(fp,file,"r");
	
	//fscanf(fp,"%*[^\n]%*c");   --need when read redmapper_dr8_public_v5.10_catalog.dat or cfht_iso_158.dat
    ngroup = 0;
	for(i=0;i<num;i++)
	{
		//fscanf(fp,"%d %lf %lf %lf",&ntmp,&tmp[0],&tmp[1],&tmp[2]);
		fscanf(fp,"%lg %lf %lf %lf",&ntmp1,&tmp[0],&tmp[1],&tmp[2]);
		pos[ngroup][0] = -tmp[0];
		pos[ngroup][1] = tmp[1];
		pos[ngroup][2] = tmp[2];
		printf("%lf %lf %f \n",pos[ngroup][0],pos[ngroup][1],pos[ngroup][2]);
		pos[ngroup][3] = DH*qromb(angu_dia_d_ww,0,tmp[2])/(1.+tmp[2]);
		printf("%lf\n",pos[ngroup][3]);
		for(j=0;j<3;j++)
			fscanf(fp,"%lf",&tmp2);
        if(pos[ngroup][2]>zbin_c[0] && pos[ngroup][2]<zbin_c[1]) ngroup++;
	}
	fclose(fp);
	printf("num=%d,ngroup=%d\n",num,ngroup);
	return ngroup;
}

SHEAR **shear_calloc(int n1,int n2)
{
	SHEAR **shear;
	int i,j;

	shear = mymalloc(sizeof(SHEAR *)*n1);
	for(i=0;i<n1;i++)
		shear[i] = mycalloc(sizeof(SHEAR),n2);
	return shear;
}

SHEARR **shearr_calloc(int n1,int n2)
{
	SHEARR **shearr;
	int i,j;

	shearr = mymalloc(sizeof(SHEARR *)*n1);
	for(i=0;i<n1;i++)
		shearr[i] = mycalloc(sizeof(SHEARR),n2);
	return shearr;
}
/*
int my_void(double (*pos)[4],int num,int zbin_use,int nxy,double zstep,int nzbin)
{
    FILE *fp;
    char file[1024];
    int i,k,tmp1=0;
    char dir[1024]="data-2017/mask_flag/maskmethod1_nlimit12_1200_rmax0.0166667";

    sprintf(file,"%s/nzbin5-ngrid1200-zstep0.1-rmax0.0166667_pos.dat",dir);
    myfopen(fp,file,"r");
    for(k=0;k<nzbin;k++)
    {
        for(i=0;i<nxy;i++)
        {
            fscanf(fp,"%lg %lg",&pos[tmp1][0],&pos[tmp1][1]);
            pos[tmp1][2] = zstep*(k+0.5);
            if(k>=zbin_use && pos[tmp1][0]!=-1000) 
            {
                pos[tmp1][3] = DH*qromb(angu_dia_d,0,pos[tmp1][2])/(1.+pos[tmp1][2]);
                printf("%lg %lg %lg %lg\n",pos[tmp1][0],pos[tmp1][1],pos[tmp1][2],pos[tmp1][3]);
                tmp1++;
            } 
            //if(k==3 && pos[tmp1][0]!=-1000) tmp1++;
        }
    }
    fclose(fp);
    return tmp1;
}
*/
int my_void(double (*pos)[5])
{
    FILE *fp;
    char file[1024];
    int i,k,nline;double tmp;
    char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/w1";
    //char dir[1024] = "data-2017/mask_flag/maskmethod1_USE";
    //sprintf(file,"%s/w1_maskmethod1_nzbin1-step0.006139-zstep0.1-zmin0.4-rmax0.025_magl-20-pos.dat",dir); 
    //sprintf(file,"%s/w1_maskmethod1_nzbin1-step0.006139-zstep0.1-zmin0.4-rmax0.0416667_magl-20-maskratio0.1-pos.dat",dir);
    sprintf(file,"%s/w1_maskmethod1_nzbin1-step0.006139-zstep0.1-zmin0.4-rmax0.0333333_magl-20.5-pos.dat",dir);
    //sprintf(file,"%s/w1_maskmethod1_nzbin1-step0.006139-zstep0.1-zmin0.4-rmax0.05_magl-20.5-maskratio0.1-pos.dat",dir);
    //sprintf(file,"%s/w1_maskmethod1_nzbin1-step0.006139-zstep0.3-zmin0.2-rmax0.0833333_magl-23-pos.dat",dir);
    //sprintf(file,"%s/w1_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.05_magl-20-pos.dat",dir);
    //sprintf(file,"%s/w1_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.0416667ANG_magl-21-pos.dat",dir);
    //sprintf(file,"%s/w1_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.0416667ANG_magl-19-maskratio0.1-pos.dat",dir);
    //sprintf(file,"%s/4field-maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.0333333_magl-18-pos.dat",dir);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg %lg",&pos[i][0],&tmp,&pos[i][1],&pos[i][2]);
        pos[i][3] = DH*qromb(angu_dia_d_ww,0,pos[i][2])/(1.+pos[i][2]);
        //printf("%lg %lg %lg %lg\n",pos[i][0],pos[i][1],pos[i][2],pos[i][3]);
    }
    fclose(fp);
    return nline;
}

int my_void_fofpos(double (*pos)[5])
{
    FILE *fp;
    char file[1024];
    int i,k,nline;double tmp;
    char dir[1024] = "data-2017/mask_flag/maskmethod1_USE/w2";
    sprintf(file,"%s/w2_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.05_magl-20-fofr0.4pos.dat",dir);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg",&pos[i][0],&pos[i][1],&pos[i][2]);
        pos[i][3] = DH*qromb(angu_dia_d_ww,0,pos[i][2])/(1.+pos[i][2]);
    }
    fclose(fp);
    return nline;
}


int random_for_void(double (*pos)[5])
{
    FILE *fp;
    int i,nline,tmp=0; 
    //char dir[1024]="data-2017/mask_flag/maskmethod1_nlimit12_1200_rmax0.0166667/w1",file[1024];
    char dir[1024]="data-2017/mask_flag/maskmethod1_USE/w1",file[1024];

    sprintf(file,"%s/random_pos_mul20.dat",dir);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg",&pos[tmp][0],&pos[tmp][1],&pos[tmp][2]);
        pos[tmp][3] = DH*qromb(angu_dia_d_ww,0,pos[tmp][2])/(1.+pos[tmp][2]);
        //printf("%lg %lg %lg %lg\n",pos[tmp][0],pos[tmp][1],pos[tmp][2],pos[tmp][3]);
        if(pos[tmp][2]>0.1) tmp++;  // z change here!!!
    }
    fclose(fp);
    return tmp;
}


int random_for_cluster(double (*pos)[5])
{
    FILE *fp;
    int i,nline,tmp=0;
    char dir[1024]="data-2017/mask_flag/maskmethod1_USE",file[1024];

    sprintf(file,"%s/4field-random-for-cluster_pos_mul10.dat",dir);
    GET_line(file,&nline);
    myfopen(fp,file,"r");
    for(i=0;i<nline;i++)
    {
        fscanf(fp,"%lg %lg %lg",&pos[tmp][0],&pos[tmp][1],&pos[tmp][2]);
        pos[tmp][3] = DH*qromb(angu_dia_d_ww,0,pos[tmp][2])/(1.+pos[tmp][2]);
        if(pos[tmp][2]>0.1) tmp++;  // z change here!!!
    }
    fclose(fp);
    return tmp;
}



/*
int read_cluster(int flagc,double (*group_pos)[5])
{
    int ngroup,ngroup_tmp=22694;//708; 
    int nxy=1200*1200,nzbin=5; double zstep=0.1;//params for void

    if(flagc == 0) ngroup_tmp=22694;
    else if(flagc == 1) ngroup_tmp=319;//26350;
    else if(flagc == -1) ngroup_tmp=nxy*nzbin;
    
    if(flagc == 0)
        ngroup = clusters(group_pos,ngroup_tmp);
    else if(flagc == 1)
        ngroup = redmapper(group_pos,ngroup_tmp);
    else if(flagc == -1)
        ngroup = my_void(group_pos);
    else if(flagc == -10)
        ngroup = my_void_fofpos(group_pos);
    else if(flagc==-11)
        ngroup = random_for_void(group_pos);
    else if(flagc == -2)
        ngroup = random_for_cluster(group_pos);
    return ngroup;
}
*/
