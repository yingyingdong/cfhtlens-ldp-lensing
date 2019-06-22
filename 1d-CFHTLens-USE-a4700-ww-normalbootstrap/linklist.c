


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "struct.h"
#include "fftw3.h"
#include "func.h"

void makell_sub(LENS *Lens,int nfit,int *hoc,int *ll,int NDIV[2],double range[2][2],double step[2],double **xgrid)
{
	int i,j,grid[2],np;
	int ndiv2=NDIV[0]*NDIV[1],index;
#define POS(i,j) Lens[i].pos[j]

	np = nfit;
	printf("creating linked list...\n");
	
	for(i=0;i<NDIV[0]*NDIV[1];i++)
		hoc[i] = -1;
	for(i=0;i<np;i++) ll[i] = -1;

	//for(i=0;i<1;i++)	
	for(i=0;i<np;i++) 
	{
		for(j=0;j<2;j++){
			//if(i<1000) printf("POS(%d,%d)=%lg,Lens[%d].pos[%d]=%f ",i,j,POS(i,j),i,j,Lens[i].pos[j]); fflush(stdout);
			grid[j] = floor((POS(i,j)-range[j][0])/step[j]);
			//if(i<1000) printf("grid[%d]=%d ",j,grid[j]);
		}
		//if(i<1000) printf("\n");
		index = grid[0]+grid[1]*NDIV[0];
		ll[i] = hoc[index];
		hoc[index] = i;
	}
/*
	int tmp = 0;
	for(i=0;i<ndiv3;i++)
		if(hoc[i] != -1) {tmp ++; index = i;}
	
	printf("tmp = %d \n",tmp);
	index = hoc[index];
	while(index != -1)
	{
		printf("Lens.z=%f ",Lens[index].pos[2]);
		index = ll[index];
	}
*/
#undef POS
	printf("finished linklist\n");
}


int *linklist_search_sub(LENS *Lens,int cenid,double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
	double range[2][2],double step[2],int NDIV[2],int nfit)
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 1;
	int index,*src;
	double cen[2],pos[2],dr; int pid;
	len = 0;
	nmax = 10000;
	for(i=0;i<2;i++) cen[i] = Lens[cenid].pos[i];
	src    = mymalloc(sizeof(int)*nmax);
	for(i=0;i<2;i++)
	{
		subbox_grid[i][0] = floor((cen[i]-rmax-range[i][0])/step[i]);
		subbox_grid[i][1] = floor((cen[i]+rmax-range[i][0])/step[i]);
		if(subbox_grid[i][0] < 0) subbox_grid[i][0] = 0;
		else if(subbox_grid[i][1] >= NDIV[i]) 
			subbox_grid[i][1] = NDIV[i]-1;
	}

	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			{
				index = i+j*NDIV[0];
				pid = hoc[index];
				//printf("index = %d ,pid = %d ,i,j,k=%d  %d  %d \n",index,pid,i,j,k);
				while(pid>=0)
				{
					for(l=0;l<2;l++) pos[l] = Lens[pid].pos[l];
					dr = distance_2dD(pos,cen);
					//printf("pid = %d,pos = %f,%f ,dr = %f \n",pid,pos[0],pos[1],dr);
					if(dr < rmax && dr >= rmin)
					{
						if(len == nmax)
						{
							nmax*=2;
							src=(int *)realloc(src,sizeof(int)*nmax);
						}
						src[len] = pid;
						len++;
					}
					pid = ll[pid];
					//if(pid >= nfit) printf("linklist wrong ! \n");
				}
				//printf("\n");
			}
			
	src = (int *)realloc(src,sizeof(int)*len);
	*Nfind = len;
	return src;
}


int *linklist_search_sub_pos(LENS *Lens,double cen[2],double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
		double range[2][2],double step[2],int NDIV[2],int nfit)
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 1;
	int index,*src;
	double pos[2],dr; int pid;
	len = 0;
	nmax = 10000;
	src    = mymalloc(sizeof(int)*nmax);
	for(i=0;i<2;i++)
	{
		subbox_grid[i][0] = floor((cen[i]-rmax-range[i][0])/step[i]);
		subbox_grid[i][1] = floor((cen[i]+rmax-range[i][0])/step[i]);
		if(subbox_grid[i][0] < 0) subbox_grid[i][0] = 0;
		else if(subbox_grid[i][1] >= NDIV[i])
			subbox_grid[i][1] = NDIV[i]-1;
	}

	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
		{
			index = i+j*NDIV[0];
			pid = hoc[index];
			while(pid>=0)
			{
				for(l=0;l<2;l++) pos[l] = Lens[pid].pos[l];
				dr = distance_2dD(pos,cen);
				if(dr < rmax && dr >= rmin)
				{
					if(len == nmax)
					{
						nmax*=2;
						src=(int *)realloc(src,sizeof(int)*nmax);
					}
					src[len] = pid;
					len++;
				}
				pid = ll[pid];
			}
		}
	src = (int *)realloc(src,sizeof(int)*len);
	*Nfind = len;
	return src;
}


int *linklist_search_sub_pos2(LENS *Lens,double cen[2],double rmax,int *Nfind,int *hoc,int *ll, \
		double range[2][2],double step[2],int NDIV[2],int nfit)
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 1;
	int index,*src;
	double pos[2],dr; int pid;
	len = 0;
	nmax = 90000;
	src    = mymalloc(sizeof(int)*nmax);
	for(i=0;i<2;i++)
	{
		subbox_grid[i][0] = floor((cen[i]-rmax-range[i][0])/step[i]);
		subbox_grid[i][1] = floor((cen[i]+rmax-range[i][0])/step[i]);
		if(subbox_grid[i][0] < 0) subbox_grid[i][0] = 0;
		else if(subbox_grid[i][1] >= NDIV[i])
			subbox_grid[i][1] = NDIV[i]-1;
	}

	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
		{
			index = i+j*NDIV[0];
			pid = hoc[index];
			while(pid>=0)
			{
				for(l=0;l<2;l++) pos[l] = Lens[pid].pos[l];
				dr = distance_2dD(pos,cen);
				if(dr < rmax)
				{
					if(len == nmax)
					{
						nmax*=2;
						src=(int *)realloc(src,sizeof(int)*nmax);
					}
					src[len] = pid;
					len++;
				}
				pid = ll[pid];
			}
		}
	src = (int *)realloc(src,sizeof(int)*len);
	*Nfind = len;
	return src;
}


int *linklist_search_sub_pos2_countzbin12(LENS *Lens,double cen[4],double rmax,int *hoc,int *ll,\
        double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,int **flag_mask,int flagd,
        double mask_ratio)
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 0,ds=0;
	int index,*src,grid,*nfind,ds0,dsi;
	double pos[2],dx,dy,dr,zmin,zmax; int pid;
    
    nfind = mycalloc(sizeof(int),nzbin);
    zmin=zmean-zstep*0.5;
    zmax=zmean+zstep*0.5;

#define POS(i,j) Lens[i].pos[j]
    for(i=0;i<2;i++)
    {
        subbox_grid[i][0] = MAX_I(0,floor((cen[i]-rmax-range[i][0])/step[i]));
        subbox_grid[i][1] = MIN_I(NDIV[i]-1,floor((cen[i]+rmax-range[i][0])/step[i]));
    }
    
    //ds0 = region(subbox_grid[0],subbox_grid[1]);
    for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
    {
        for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
        {
            //if(flag_mask[j][i]==0) {flag=1; break;}
            ds++; if(flag_mask[j][i]==0) {flag ++;  continue;}
            
            index = i+j*NDIV[0];
            pid = hoc[index];
            while(pid>=0)
            {
                if(POS(pid,2) < zmax && POS(pid,2) > zmin)
                {
                    if(flagd==0)
                    {
                        dx = POS(pid,0) - cen[0];
                        dy = POS(pid,1) - cen[1];
                    }
                    dr = sqrt(dx*dx + dy*dy);
                    if(dr < rmax) 
                        nfind[(int)floor( (POS(pid,2)-zmin)/zstep )] ++;
                }
                pid = ll[pid];
            }
        }
        //if(flag == 1) break;
    }
    //if(flag > 0) memset(nfind,-1,sizeof(int)*nzbin);
    if(flag*1./ds > mask_ratio) memset(nfind,-1,sizeof(int)*nzbin);
    return nfind;
#undef POS
}

int gain_rbin(int flag_lg,double dr,double rmin,double rstep)
{
    if(flag_lg == 0) return floor(dr/rstep);
    else return floor( (log10(dr)-log10(rmin))/rstep );
}

int ***linklist_search_sub_pos2_countzbin12_record(LENS *Lens,double cen[5],double rmax,double rmin,double rstep,
        int *hoc,int *ll,double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,
        int **flag_mask,int flagd,int flag_lg,double mask_ratio,int *flag_use,int ntr[2])
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 0,ds=0;
	int index,*src,grid,ds0,dsi,rtmp,***nfind_tmp;
	double pos[2],dx,dy,dr,zmax,zmin; int pid;
   
    nfind_tmp = I_malloc3(nzbin,ntr[0],ntr[1]);
    for(i=0;i<nzbin;i++) for(j=0;j<ntr[0];j++) for(k=0;k<ntr[1];k++) 
        nfind_tmp[i][j][k]=0; //zero the arr
    zmin=zmean-zstep*0.5;
    zmax=zmean+zstep*0.5;

#define POS(i,j) Lens[i].pos[j]
    for(i=0;i<2;i++)
    {
        subbox_grid[i][0] = MAX_I(0,floor((cen[i]-rmax-range[i][0])/step[i]));
        subbox_grid[i][1] = MIN_I(NDIV[i]-1,floor((cen[i]+rmax-range[i][0])/step[i]));
    }
    
    //ds0 = region(subbox_grid[0],subbox_grid[1]);
    for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
    {
        for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
        {
            //if(flag_mask[j][i]==0) {flag=1; break;}
            ds++; if(flag_mask[j][i]==0) {flag ++;  continue;}  //change here!!!
            
            index = i+j*NDIV[0];
            pid = hoc[index];
            while(pid>=0)
            {
                if(POS(pid,2) < zmax && POS(pid,2) > zmin)
                {
                    dx = POS(pid,0) - cen[0];
                    dy = POS(pid,1) - cen[1];
                    dr = sqrt(dx*dx + dy*dy);
                    if(dr < rmax && dr >= rmin)
                    {
                        if(flagd == 1) dr *=(cen[3]*au*(1+cen[2])); //comoving distance
                        rtmp = gain_rbin(flag_lg,dr,rmin,rstep);
                        if(rtmp < ntr[1])
                        {
                            nfind_tmp[(int)floor((POS(pid,2)-zmin)/zstep)][0][rtmp] ++;
                        }
                    }
                }
                pid = ll[pid];
            }
        }
        //if(flag == 1) break;
    }
    
    if(flag*1./ds <= mask_ratio) *flag_use = 1;
    else *flag_use = 0;
    return nfind_tmp;
#undef POS
}

int ***linklist_search_sub_pos2_countzbin12_record_nomask(LENS *Lens,double cen[5],double rmax,double rmin,double rstep,
        int *hoc,int *ll,double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,
        int flagd,int flag_lg,int *flag_use,int ntr[2])
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 0,ds=0;
	int index,*src,grid,ds0,dsi,rtmp,***nfind_tmp;
	double pos[2],dx,dy,dr,zmax,zmin; int pid;
   
    nfind_tmp = I_malloc3(nzbin,ntr[0],ntr[1]);
    for(i=0;i<nzbin;i++) for(j=0;j<ntr[0];j++) for(k=0;k<ntr[1];k++) 
        nfind_tmp[i][j][k]=0; //zero the arr
    zmin=zmean-zstep*0.5;
    zmax=zmean+zstep*0.5;

#define POS(i,j) Lens[i].pos[j]
    for(i=0;i<2;i++)
    {
        subbox_grid[i][0] = MAX_I(0,floor((cen[i]-rmax-range[i][0])/step[i]));
        subbox_grid[i][1] = MIN_I(NDIV[i]-1,floor((cen[i]+rmax-range[i][0])/step[i]));
    }
    
    //ds0 = region(subbox_grid[0],subbox_grid[1]);
    for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
    {
        for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
        {
            
            index = i+j*NDIV[0];
            pid = hoc[index];
            while(pid>=0)
            {
                if(POS(pid,2) < zmax && POS(pid,2) > zmin)
                {
                    dx = POS(pid,0) - cen[0];
                    dy = POS(pid,1) - cen[1];
                    dr = sqrt(dx*dx + dy*dy);
                    if(dr < rmax && dr >= rmin)
                    {
                        if(flagd == 1) dr *=(cen[3]*au*(1+cen[2])); //comoving distance
                        rtmp = gain_rbin(flag_lg,dr,rmin,rstep);
                        if(rtmp < ntr[1])
                        {
                            nfind_tmp[(int)floor((POS(pid,2)-zmin)/zstep)][0][rtmp] ++;
                        }
                    }
                }
                pid = ll[pid];
            }
        }
        //if(flag == 1) break;
    }
    
    return nfind_tmp;
#undef POS
}

int *linklist_search_sub_pos2_countzbin12_COM(LENS *Lens,double cen[4],double rmax,int *hoc,int *ll,\
        double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,int **flag_mask,
        double mask_ratio,double DCOM[])
{
	int i,j,k,l,subbox_grid[2][2],nmax,len,flag = 0,ds=0;
	int index,*src,grid,*nfind,ds0,dsi,nzbini,maskid[2];
	double pos[2],dx,dy,dr,zmax,zmin; int pid;
    
    nfind = mycalloc(sizeof(int),nzbin);
    zmin=zmean-zstep*0.5;
    zmax=zmean+zstep*0.5;

#define POS(i,j) Lens[i].pos[j]
    for(i=0;i<2;i++)
    {
        subbox_grid[i][0] = MAX_I(0,floor((cen[i]-rmax-range[i][0])/step[i]));
        subbox_grid[i][1] = MIN_I(NDIV[i]-1,floor((cen[i]+rmax-range[i][0])/step[i]));
    }
    
    //ds0 = region(subbox_grid[0],subbox_grid[1]);
    for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
    {
        for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
        {
            //if(flag_mask[j][i]==0) {flag=1; break;}
            ds++; if(flag_mask[j][i]==0) {flag ++;  continue;}
            
            index = i+j*NDIV[0];
            pid = hoc[index];
            while(pid>=0)
            {
                if(POS(pid,2) < zmax && POS(pid,2) > zmin)
                {
                    dx = POS(pid,0) - cen[0];
                    dy = POS(pid,1) - cen[1];
                    dr = sqrt(dx*dx + dy*dy); 
                    nzbini = floor( (POS(pid,2)-zmin)/zstep );
                    if(dr < rmax*DCOM[0]/DCOM[nzbini]) 
                        nfind[nzbini] ++;
                }
                pid = ll[pid];
            }
        }
        //if(flag == 1) break;
    }
    if(flag > 0) memset(nfind,-1,sizeof(int)*nzbin);
    //if(flag*1./ds > mask_ratio) memset(nfind,-1,sizeof(int)*nzbin);
    return nfind;
#undef POS
}


void check_num(LENS *Lens,int NDIV[2],double **xgrid,int *hoc,int *ll)
{
   	int ndiv2=NDIV[0]*NDIV[1],ndiv3=NDIV[0]*NDIV[1]*NDIV[2],j,k,l;
	int tmp = 0,tmp1 = 0,i,pid,index,use[3];
	
	for(i=0;i<NDIV[0];i++)
		for(j=0;j<NDIV[1];j++)
			{
				index = i+j*NDIV[0];
				pid = hoc[index];
				//while(pid >= 0 && tmp < 20)
				while(pid >= 0)
				{
					tmp ++;
					use[0] = i;
					use[1] = j;
					for(l=0;l<2;l++) if(Lens[pid].pos[l] < xgrid[l][use[l]] || 
							Lens[pid].pos[l] > xgrid[l][use[l]+1])
					{
						printf("grid wrong !,xgrid[%d][%d]=%f,Lens[%d].pos[%d]=%f \n",l,use[l],xgrid[l][use[l]],pid,l,Lens[pid].pos[l]);
						tmp1 ++;
					}
					if(Lens[pid].e1 == 0 && Lens[pid].e2 == 0) printf("wrong ! wt=%f \n",Lens[pid].wt);
					pid = ll[pid];
				}
			}
	printf("tmp = %d\n",tmp);
}


void check_search_sub(int Nfound,int *found,int nfit)
{
	int i;
	
	for(i=0;i<Nfound;i++) 
		if(found[i] >= nfit) printf("check in search sub is wrong ! \n");
}

/*
int *linklist_search_sub(SUBCATALOGUE SubTb,int cenid,float rmax,float (*subpos)[3],int *Nfind,int *hoc,int *ll)
{
	int i,j,k,subbox_grid[3][2],nmax,len;
	int ndiv2 = NDIV*NDIV,index,*src;
	HBTReal *cen,dr; int pid;

	len = 0;
	nmax = 10;
	cen    = subpos[cenid];
	src    = mymalloc(sizeof(int)*nmax);
	for(i=0;i<3;i++)
	{
		subbox_grid[i][0] = floor((cen[i]-rmax-range[i][0])/step[i]);
		subbox_grid[i][1] = floor((cen[i]+rmax-range[i][0])/step[i]);
	}

	for(i=subbox_grid[0][0];i<subbox_grid[0][1]+1;i++)
		for(j=subbox_grid[1][0];j<subbox_grid[1][1]+1;j++)
			for(k=subbox_grid[2][0];k<subbox_grid[2][1]+1;k++)
			{
				index = ijk_shift(i)+ijk_shift(j)*NDIV+ijk_shift(k)*ndiv2;
				pid = hoc[index];
				while(pid>0)  //? >=0
				{
					dr = distance(subpos[pid],cen);
					if(pid!=cenid && dr<=rmax && SubTb.SubLen[pid]>0)
					{
						if(len == nmax)
						{
							nmax*=2;
							src=(int *)realloc(src,sizeof(int)*nmax);
						}
						src[len] = pid;
						len++;
					}

					pid = ll[pid];
				}
			}
	src = (int *)realloc(src,sizeof(int)*len);
	*Nfind = len;
	return src;
}
*/

