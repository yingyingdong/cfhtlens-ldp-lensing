#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "struct_jz.h"
#include "fftw3.h"
#include "func.h"
#include "fft_m.h"


fftw_complex *gain_ll_num(int nx,int ny,int *hoc,int *ll)
{
    int i,j,index,pid;
    fftw_complex *ima;
    ima = zero_ima(nx,ny);

    for(i=0;i<nx;i++) //lie shu, x axis
        for(j=0;j<ny;j++) //hang shu, y axis
        {
            index = i+j*nx;
            pid = hoc[index];
            while(pid>=0)
            {
                ima[ELEM(j,i,nx)][0] +=1;
                pid = ll[pid];
            }
        }
    return ima;
}

fftw_complex *make_convolve_fft(fftw_complex *ima,fftw_complex *w,int nx,int ny,fftw_plan p)
{
    fftw_complex *w_f,*ima_f,*c_f,*c_ff;

    w_f = make_2dfft(w,nx,ny,p,0);
    ima_f = make_2dfft(ima,nx,ny,p,0);
    c_f = fft_convolve(ima_f,w_f,nx,ny);
    c_ff = make_2difft(c_f,nx,ny,p);
    normal_fftifft(c_ff,nx,ny);
    return fftshift_2d(c_ff,nx,ny);
}

double gauss_here(double x,double y,double a)
{
    return 1./(2*pi*a*a)*exp( -(x*x+y*y)/(2*a*a) );
}

fftw_complex *gain_gauss(double a,double **xgrid,int nx,int ny,double xc,double yc)
{
    int i,j;
    fftw_complex *w;
    w = zero_ima(nx,ny);
    for(i=0;i<ny;i++) //hang ,y axis
        for(j=0;j<nx;j++) //lie ,x axis
        {
             w[ELEM(i,j,nx)][0] = gauss_here(xgrid[0][j]-xc , xgrid[1][i]-yc,a);
        }
    return w;
}

void gain_name_here(char *nameima,char *nameconv,char *nameimaconv,int nfield)
{
    sprintf(nameima,"data-2017/mask_flag/maskmethod1_USE/w%d/ima.dat",nfield+1);
    sprintf(nameconv,"data-2017/mask_flag/maskmethod1_USE/w%d/conv.dat",nfield+1);
    sprintf(nameimaconv,"data-2017/mask_flag/maskmethod1_USE/w%d/ima_conv.dat",nfield+1);

}
int **mask_region(int nx,int ny,int *hoc,int *ll,double step[2],
        double **xgrid,double xc,double yc,int mask_method,int nfield)
{
    int i,j,**flag; double ds,a;
    char nameima[512],nameconv[512],nameimaconv[512];
    fftw_complex *w,*ima,*ima_c;
    fftw_plan p;

    gain_name_here(nameima,nameconv,nameimaconv,nfield);
    flag = I_ones2(ny,nx);
    ima = gain_ll_num(nx,ny,hoc,ll);
    a = MIN_D(step[0],step[1]); printf("step[0]=%f,step[1]=%f,a=%e\n",step[0],step[1],a);
    w = gain_gauss(a,xgrid,nx,ny,xc,yc);
    normal_ima(w,nx,ny);   //normallize the window function 
    ds = step[0]*step[1]*60*60;
    
    ima_c = make_convolve_fft(ima,w,nx,ny,p);
    save_complex2d(ima,nx,ny,nameima);
    save_complex2d(w,nx,ny,nameconv);
    save_complex2d(ima_c,nx,ny,nameimaconv);
    for(i=0;i<ny;i++) //hang, y axis
        for(j=0;j<nx;j++) //lie, x axis
        {
            if(mask_method == 1 && ima_c[ELEM(i,j,nx)][0]/ds < 12)
                flag[i][j] = 0;
            else if(mask_method == 0 && ima[ELEM(i,j,nx)][0] == 0)
                flag[i][j] = 0;
        }
    return flag;
}

int **load_flag_mask(int nfield,int nx,int ny,int load_mask)
{
    int i,j,**flag;
    char file[1024];
    
    if(load_mask==1)
        sprintf(file,"data-2017/mask_flag/maskmethod1_USE/w%d/w%d_maskmethod1_step0.006139_voidflag.dat",nfield+1,nfield+1);
    if(load_mask==2)
        sprintf(file,"data-2017/mask_flag/2load_cfht/w%d_xy_mask.dat",nfield+1);
    flag = load_int2d(nx,ny,file);
    return flag;
}


