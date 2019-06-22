#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "fitsio.h"
#include "fftw3.h"
#include "func.h"
#include "fft_m.h"
#include "omp.h"

int ELEM(int r,int c,int N)
{
    return r*N+c;
}

fftw_complex *zero_ima(int nx,int ny)
{
    fftw_complex *ima;
    ima = mymalloc_C(sizeof(fftw_complex)*nx*ny);
    return memset(ima,0,sizeof(fftw_complex)*nx*ny);
}

float sum_2d(fftw_complex *arr,int nx,int ny)
{
    int i; float sum=0;
    for(i=0;i<nx*ny;i++) sum +=arr[i][0];
    return sum;
}

void swap(fftw_complex *v1, fftw_complex *v2)
{
    fftw_complex tmp;

    tmp[0] = *v1[0]; tmp[1] = *v1[1];
    *v1[0] = *v2[0]; *v1[1] = *v2[1];
    *v2[0] = tmp[0]; *v2[1] = tmp[1];
}

double max_2d_complex(fftw_complex *arr,int nx,int ny,int flag)
{
    int i,j;
    double max;

    max =arr[ELEM(0,0,ny)][flag];
    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
        {
           if(arr[ELEM(i,j,nx)][flag] > max)
               max = arr[ELEM(i,j,nx)][flag];
        }
    return(max);
}


double pow_complex(fftw_complex image)
{
    return image[0]*image[0]+image[1]*image[1];
}

void shift_phase_2d(fftw_complex *data, int nx,int ny)
{
    int i,j;

    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
        {
            data[ELEM(i,j,nx)][0] *=pow(-1.,i+j);
            data[ELEM(i,j,nx)][1] = 0;
        }
}

fftw_complex *make_2dfft(fftw_complex *ima,int nx,int ny,fftw_plan p,int phase)
{
    fftw_complex *ima_f;
    
    if(phase == 1) shift_phase_2d(ima,nx,ny);
    ima_f = mymalloc_C(sizeof(fftw_complex)*nx*ny);
    p=fftw_plan_dft_2d(ny,nx,ima,ima_f,\
            FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    return ima_f;
}

fftw_complex *make_2difft(fftw_complex *ima,int nx,int ny,fftw_plan p)
{
    fftw_complex *ima_f;

    ima_f = mymalloc_C(sizeof(fftw_complex)*nx*ny);
    p=fftw_plan_dft_2d(ny,nx,ima,ima_f,\
            FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    return ima_f;
}

void normal_fftifft(fftw_complex *ima,int nx,int ny)
{
    int i,j;

    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
        {
            ima[ELEM(i,j,nx)][0] /=(nx*ny)*1.;
        }
}


void normal_ima(fftw_complex *ima,int nx,int ny)
{
    int i; float sum_ima;
    sum_ima = sum_2d(ima,nx,ny);
    for(i=0;i<nx*ny;i++)
        ima[i][0] /= sum_ima*1.;
}

void pow_2darr(fftw_complex *ima_f,int nx,int ny)
{
    int i,j;

    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
            ima_f[ELEM(i,j,nx)][0] = pow_complex(ima_f[ELEM(i,j,nx)]);
}

double real_multi_fft(fftw_complex arr1,fftw_complex arr2)
{
    return arr1[0]*arr2[0]-arr1[1]*arr2[1];
}

double imag_multi_fft(fftw_complex arr1,fftw_complex arr2)
{
    return arr1[0]*arr2[1]+arr1[1]*arr2[0];
}


fftw_complex *fft_convolve(fftw_complex *image_f,fftw_complex *psf_f,int nx,int ny)
{
    int i,j;
    fftw_complex *c_f;
    c_f = mymalloc_C(sizeof(fftw_complex)*nx*nx);

    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++)  //lie
        {
            c_f[ELEM(i,j,nx)][0] = real_multi_fft(image_f[ELEM(i,j,nx)],psf_f[ELEM(i,j,nx)]);
            c_f[ELEM(i,j,nx)][1] = imag_multi_fft(image_f[ELEM(i,j,nx)],psf_f[ELEM(i,j,nx)]);
        }
    return c_f;
}

fftw_complex *deconvolve(fftw_complex *image_f,fftw_complex *psf_f,fftw_complex *gauss_f,int nx)
{
    int i,j; double max=0;;
    fftw_complex *image_t_f;
    image_t_f = mymalloc_C(sizeof(fftw_complex)*nx*nx);
    for(i=0;i<nx*nx;i++)
    {
        image_t_f[i][0] =
            sqrt(pow_complex(image_f[i])*(pow_complex(gauss_f[i])/pow_complex(psf_f[i])));
        image_t_f[i][1] = 0;
    }
    return image_t_f;
}

double Max_fftw_2darr(fftw_complex *ima,int nx,int ny)
{
    int i,j;
    double max=0.;
    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
            if(ima[ELEM(i,j,nx)][0] > max) 
                max = ima[ELEM(i,j,nx)][0];    
    return max;
}

double area_thresh(fftw_complex *ima,double thresh,int nx,int ny)
{
    int i,j;
    double s=0;
    
    for(i=0;i<ny;i++) //hang
        for(j=0;j<nx;j++) //lie
            if(ima[ELEM(i,j,nx)][0] > thresh) s+=1;
    return s;
}

double poisson_thresh(fftw_complex *ima,double thresh,int nx,int ny)
{
    int i,j;
    double s=0.,flux=0.;

    for(i=0;i<ny;i++)  //hang
        for(j=0;j<nx;j++) //lie
            if(ima[ELEM(i,j,nx)][0] > thresh){
                s +=1;
                flux +=ima[ELEM(i,j,nx)][0];
            }
    return flux/sqrt(s);
}

void gain_gauss_fft(fftw_complex *psf_f,fftw_complex *gal_f,fftw_complex *noise_f,
        int nx,float *kx,float *gamma1,float *gamma2,float *gamma3,float h[2])
{
    int i,j;
    double thresh1,thresh2,max,s,g1=0,g2=0,de=0,h1=0,h2=0;
    double kxi,kyj,k2,ks,ks_2,ff,temp,temp1;

    max = Max_fftw_2darr(psf_f,nx,nx);
    thresh1 = max*exp(-1.);
    s = area_thresh(psf_f,thresh1,nx,nx);
    ks = sqrt(s/pi);
    ks_2 = pow(ks*psf_ratio,-2.);

    thresh2 = max*1e-5;
    //for(i=0;i<nx;i++)
    for(j=0;j<nx;j++)
    {
        kxi = kx[j];
        //for(j=0;j<nx;j++)
        for(i=0;i<nx;i++)
        {
            kyj = kx[i];
            k2 = (kxi*kxi+kyj*kyj);
            if(psf_f[ELEM(i,j,nx)][0] > thresh2)
            {
                ff = k2*ks_2;
                temp = exp(-ff)/psf_f[ELEM(i,j,nx)][0];
                temp1 = temp*(gal_f[ELEM(i,j,nx)][0]-noise_f[ELEM(i,j,nx)][0]);
                //temp1 = temp*(gal_f[ELEM(i,j,nx)][0]);
                g1 -=temp1*(kxi*kxi-kyj*kyj);
                g2 -=temp1*2*kxi*kyj;
                de +=temp1*k2*(2.-ff);
                h1 +=temp1*ks_2*(k2*k2-8.*kxi*kxi*kyj*kyj);
                h2 +=temp1*ks_2*4.*kxi*kyj*(kxi*kxi-kyj*kyj);
            }
        }
    }
    *gamma1 = g1;
    *gamma2 = g2;
    *gamma3 = de;
    h[0] = h1;   h[1] = h2;
}



fftw_complex *fftshift_2d(fftw_complex *data, int nx,int ny)
{
    int i,j = 0;
    int cx = nx/2,cy=ny/2;
    fftw_complex *c_s,*c_ss,tmp;
    c_s = mymalloc_C(sizeof(fftw_complex)*nx*ny);
    c_ss = mymalloc_C(sizeof(fftw_complex)*nx*ny);

    //printf("ncount=%d\n",c);
    // For odd and for even numbers of element use different algorithm
    for(i=0;i<ny;i++)  //hang
        for(j=0;j<cx;j++) //lie
        {
            c_s[ELEM(i,j,nx)][0] = data[ELEM(i,j+cx,nx)][0];
            c_s[ELEM(i,j,nx)][1] = data[ELEM(i,j+cx,nx)][1];
            c_s[ELEM(i,j+cx,nx)][0] = data[ELEM(i,j,nx)][0];
            c_s[ELEM(i,j+cx,nx)][1] = data[ELEM(i,j,nx)][1];
        }
    
    for(i=0;i<cy;i++) //hang
        for(j=0;j<nx;j++) //lie
        {
            c_ss[ELEM(i,j,nx)][0] = c_s[ELEM(i+cy,j,nx)][0];
            c_ss[ELEM(i,j,nx)][1] = c_s[ELEM(i+cy,j,nx)][1];
            c_ss[ELEM(i+cy,j,nx)][0] = c_s[ELEM(i,j,nx)][0];
            c_ss[ELEM(i+cy,j,nx)][1] = c_s[ELEM(i,j,nx)][1];
        }
    return c_ss;
}



/*
void ifftshift(fftw_complex *data, int count)
{
    int k = 0;
    int c = (int) floor((float)count/2);
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
            swap(&data[k], &data[k+c]);
    }
    //else
    //{
    //    fftw_complex tmp = data[count - 1];
    //    for (k = c-1; k >= 0; k--)
    //    {
    //        data[c + k + 1] = data[k];
    //        data[k] = data[c + k];
    //    }
    //    data[c] = tmp;
    //}
}
*/

float Moffat(float rij,float rd,float rc)
{
    return pow((1.+(rij/rd)*(rij/rd)),-3.5);
}

void Moffat_profile(fftw_complex *psf,int npsf,int nx,float rd,float rc,float ncx,float ncy)
{
    int i,j;
    float r;
    for(i=0;i<nx;i++)
        for(j=0;j<nx;j++)
        {
            r = sqrt((i-ncx+0.5)*(i-ncx+0.5)+(j-ncy+0.5)*(j-ncy+0.5));
            if(r < rc) psf[ELEM(i,j,nx)][0] = Moffat(r,rd,rc);  //not rc?
            else psf[ELEM(i,j,nx)][0] = 0;
        }
}


void write_imagepow_txt(fftw_complex *image,int nx,int ny,char *name)
{
    int i,j,k;
    FILE *fp;

    myfopen(fp,name,"w");
    for(i=0;i<ny;i++) //hang
    {
        for(j=0;j<nx;j++) //lie
            fprintf(fp,"%f ",sqrt(pow_complex(image[ELEM(i,j,nx)])));
        fprintf(fp,"\n");
    }
    fclose(fp);
}


/*
void add_gauss_noise(struct drand48_data *drand_buf,fftw_complex *image,
        fftw_complex *noise,int nx,float sn)
{
    int i,j;
    double mean=0,max,thresh,poisson,sigma;

    max = max_2d_complex(image,nx,nx,0);
    thresh = max*exp(-1.);
    poisson = poisson_thresh(image,thresh,nx,nx);
    for(i=0;i<nx;i++)
        for(j=0;j<nx;j++)
        {
            sigma = poisson/sn;
            //noise = mean + sigma*gaussian_std_normal(drand_buf);
            image[ELEM(i,j,nx)][0] +=mean + sigma*gaussian_std_normal(drand_buf);
            noise[ELEM(i,j,nx)][0]  =mean + sigma*gaussian_std_normal(drand_buf);
        }
}
*/

void save_complex2d(fftw_complex *arr,int nx,int ny,char *name)
{
    FILE *fp;
    int i,j;
    myfopen(fp,name,"w");

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fprintf(fp,"%e ",arr[ELEM(i,j,nx)][0]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}
