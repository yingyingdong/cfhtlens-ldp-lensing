#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fftw3.h"
#include "func.h"


double MAX_D(double a,double b)
{
    if(a>b) return a;
    else return b;
}

float MAX_F(float a,float b)
{
	if(a>b) return a;
	else return b;
}

int MAX_I(int a,int b)
{
	if(a>b) return a;
	else return b;
}

int MIN_I(int a,int b)
{
	if(a<b) return a;
	else return b;
}

double MIN_D(double a,double b)
{
    if(a<b) return a;
    else return b;
}

void *mymalloc(size_t n)
{void * mem;
	if(n)
	{
		if(!(mem=malloc(n)))
		{printf("failed to allocate memory for %u bytes.\n",(unsigned) n);fflush(stdout);
			exit(1);
		}
	}
	else
	{
		mem=NULL;
	}
	return mem;
}

void *mycalloc(size_t n,int i)
{void * mem;
	if(n)
	{
		if(!(mem=calloc(i,n)))
		{printf("failed to allocate memory for %u bytes.\n",(unsigned) n);fflush(stdout);
			exit(1);
		}
	}
	else
	{
		mem=NULL;
	}
	return mem;
}


void *mymalloc_C(size_t n)  //complex in fftw3
{
	void * mem;
	if(n)
	{
		if(!(mem=fftw_malloc(n)))
		{printf("failed to allocate memory for %u bytes.\n",(unsigned) n);fflush(stdout);
			exit(1);
		}
	}
	else
	{
		mem=NULL;
	}
	return mem;
}

int **I_malloc2(int n1,int n2)   //flag=0: mymalloc  ,flag=1: fftw_malloc
{   
	int i;
	int **arr;

	arr=mymalloc(sizeof(int *)*n1);
	for(i=0;i<n1;i++) arr[i]=mymalloc(sizeof(int)*n2);
	return arr;
}   

int *I_ones(int n1)
{
    int i;
    int *arr;
    arr=mymalloc(sizeof(int)*n1);
    for(i=0;i<n1;i++) arr[i] = 1;
    return arr;
}

int **I_ones2(int n1,int n2)   //flag=0: mymalloc  ,flag=1: fftw_malloc
{
    int i,j;
    int **arr;

    arr=mymalloc(sizeof(int *)*n1);
    for(i=0;i<n1;i++) //hang
    {
        arr[i] = mymalloc(sizeof(int)*n2);
        for(j=0;j<n2;j++) //lie
            arr[i][j] = 1;
    }
    return arr;
}

int ***I_malloc3(int n1,int n2,int n3)   //flag=0: mymalloc  ,flag=1: fftw_malloc
{
    int i,j;
    int ***arr;

    arr=mymalloc(sizeof(int **)*n1);
    for(i=0;i<n1;i++)
        arr[i] = I_malloc2(n2,n3);
    return arr;
}


float **F_malloc2(int n1,int n2)
{
	int i;
	float **arr;

	arr=mymalloc(sizeof(float *)*n1);
	for(i=0;i<n1;i++) arr[i]=mymalloc(sizeof(float)*n2);
	return arr;
}

int **I_calloc2(int n1,int n2)
{   
	int i;
	int **arr;

	arr=mymalloc(sizeof(int *)*n1);
	for(i=0;i<n1;i++) arr[i]=mycalloc(sizeof(int),n2);
	return arr;
}

int ***I_calloc3(int n1,int n2,int n3)
{
    int i,***arr;

    arr=mymalloc(sizeof(int **)*n1);
    for(i=0;i<n1;i++) arr[i]=I_calloc2(n2,n3);
    return arr;
}


float **F_calloc2(int n1,int n2)
{
	int i;
	float **arr;

	arr=mymalloc(sizeof(float *)*n1);
	for(i=0;i<n1;i++) arr[i]=mycalloc(sizeof(float),n2);
	return arr;
}

double **D_calloc2(int n1,int n2)
{
	int i;
	double **arr;

	arr=mymalloc(sizeof(double *)*n1);
	for(i=0;i<n1;i++) arr[i]=mycalloc(sizeof(double),n2);
	return arr;
}

double **D_malloc2(int n1,int n2)
{
	int i;
	double **arr;

	arr=mymalloc(sizeof(double *)*n1);
	for(i=0;i<n1;i++) arr[i]=mymalloc(sizeof(double)*n2);
	return arr;
}


fftw_complex **C_malloc2(int n1,int n2)
{
	int i;
	fftw_complex **arr;
	
	arr=(fftw_complex **) fftw_malloc(sizeof(fftw_complex *)*n1);
	for(i=0;i<n1;i++) arr[i]=\
		(fftw_complex *) fftw_malloc(sizeof(fftw_complex)*n2);
	return arr;	
}

void myfree(void *mem)
{
	if(mem!=NULL)
		free(mem);
}

void free_I_malloc2(int **arr,int n)
{
	int i;

	for(i=0;i<n;i++) free(((int **)arr)[i]);
	free(arr);
}

void free_D_malloc2(double **arr,int n)
{
	int i;

	for(i=0;i<n;i++) free(((double **)arr)[i]);
	free(arr);
}


void free_F_malloc2(float **arr,int n)
{
	int i;

	for(i=0;i<n;i++) free(((float **)arr)[i]);
	free(arr);
}
/*
void free_C_malloc2(fftw_complex **arr,int n)
{
	int i;

	for(i=0;i<n;i++) fftw_free(((fftw_complex **)arr)[i]);
	fftw_free(arr);
}
*/


double ***D_malloc3(int n1,int n2,int n3)
{
	int i,j,k;
	double ***arr;

	arr = mymalloc(sizeof(double **)*n1);
	for(i=0;i<n1;i++)
		arr[i] = mymalloc(sizeof(double *)*n2);
	for(i=0;i<n1;i++)
		for(j=0;j<n2;j++)
			arr[i][j] = mymalloc(sizeof(double)*n3);

	return arr;
}




int sum_2d_real_I(int **arr,int nx,int ny)
{
	int i,j;
	int sum=0;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			sum +=arr[i][j];
	return(sum);
}

float sum_2d_real_F(float **arr,int nx,int ny)
{
	int i,j;
	double sum=0;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			sum +=arr[i][j];
	return(sum);
}


double sum_2d_real_D(double **arr,int nx,int ny)
{
	int i,j;
	double sum=0;
	
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			sum +=arr[i][j];
	return(sum);
}

double mean_2d_real_D(double **arr,int nx,int ny)
{
	int i,j;
	double mean=0;

	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
			mean +=arr[i][j]/(nx*ny);
	return(mean);
}

double sum_1d_real(void *arr,int n,int flag)
{
	int i;
	double sum=0;

	for(i=0;i<n;i++)
	{
		if(flag==0) sum +=((int *)arr)[i];
		if(flag==1) sum +=((float *)arr)[i];
	}
	return(sum);
}
/*
void initial_1d(void *arr,int n,int flag)
{
	int i;

	for(i=0;i<n;i++){
		if(flag==0) ((int *)arr)[i]=0;
		if(flag==1) ((float *)arr)[i]=0;
		if(flag==2) {
			((fftw_complex *)arr)[i][0]=0;
			((fftw_complex *)arr)[i][1]=0;
		}
	}
}
*/
int *sort_shift(int nsize)
{
	int i,*arr,N21=nsize/2+1,*tmp;

	printf("N21=%d \n",N21);
	arr=mymalloc(sizeof(int)*nsize);
	for(i=0;i<nsize;i++)
	{
		if(i<N21-1 && i+N21<=nsize-1) arr[i]=i+N21;
		else arr[i]=i-(nsize-(N21));
	//	printf("arr[%d]=%d ",i,arr[i]);
	}
	return arr;
}

void shift_1d_F(float *arr,int nsize,int *sort)
{
	int i;
	float *tmp;

	tmp = mymalloc(sizeof(float)*nsize);
	for(i=0;i<nsize;i++)
		tmp[i]=arr[i];
	
	for(i=0;i<nsize;i++)
	{	arr[i]=tmp[sort[i]];	
//		printf("index=%d ,%f \n",sort[i],arr[i]);
	}
}
/*
void shift_1d_C(fftw_complex *arr,int nsize,int *sort)
{
	int i;
	float (*tmp)[2];

	tmp = mymalloc(sizeof(fftw_complex)*nsize*2);
	for(i=0;i<nsize;i++){
		tmp[i][0]=arr[i][0];
		tmp[i][1]=arr[i][1];
	}

	for(i=0;i<nsize;i++)
	{
		arr[i][0]=tmp[sort[i]][0];
		arr[i][1]=tmp[sort[i]][1];
	}
}

*/
float rotation_x(float theta,float x,float y)
{
	return cos(theta)*x-sin(theta)*y;
}

float rotation_y(float theta,float x,float y)
{
	return sin(theta)*x+cos(theta)*y;
}

float distance_2dF(float x[2],float y[2])
{
    float dx[2];
    dx[0]=x[0]-y[0];
    dx[1]=x[1]-y[1];  
    return sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
}

double distance_2dD(double x[2],double y[2])
{   
	double dx[2];
	dx[0]=x[0]-y[0];
	dx[1]=x[1]-y[1];  
	return sqrt(dx[0]*dx[0]+dx[1]*dx[1]);
}   

void save_int2d(int **arr,int nx,int ny,char *name)
{
    FILE *fp;
    int i,j;
    myfopen(fp,name,"w");

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fprintf(fp,"%d ",arr[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
}

int **load_int2d(int nx,int ny,char *name)
{
    FILE *fp;
    int i,j,**arr;
    arr = I_malloc2(ny,nx);
    myfopen(fp,name,"r");

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fscanf(fp,"%d",&arr[i][j]);
    }
    fclose(fp);
    return arr;
}

void save_double2d(double **arr,int nx,int ny,char *name)
{
    FILE *fp;
    int i,j;
    myfopen(fp,name,"w");

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fprintf(fp,"%lg ",arr[i][j]); //err %e
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void save_double2d4(double (*arr)[4],int nx,int ny,char *name)
{
    FILE *fp;
    int i,j;
    myfopen(fp,name,"w");

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fprintf(fp,"%lg ",arr[i][j]); //err %e
        fprintf(fp,"\n");
    }
    fclose(fp);
}

double **load_double2d(int nx,int ny,char *name)
{
    FILE *fp;
    int i,j;double **arr;
    myfopen(fp,name,"r");
    arr = D_malloc2(ny,nx);

    for(i=0;i<ny;i++)
    {
        for(j=0;j<nx;j++)
            fscanf(fp,"%lf",&arr[i][j]);
    }
    fclose(fp);
    return arr;
}

void GET_line(char *file,int *nobj)
{
    FILE *fp; int n = 0;
    char buf[1024];
    myfopen(fp,file,"r");

    while(fgets(buf,1024,fp) != NULL)
    {
        if(buf[strlen(buf) - 1] == '\n')
            n++;
    }
    fclose(fp);
    *nobj = n;
    printf("nobj = %d \n",*nobj);

}

double *return_drand48(int nrand)
{
    int i,seed;
    static struct drand48_data drand_buf;
    double *random;
    static int iset_drand48=0;

    random = mymalloc(sizeof(double)*nrand);
    if(iset_drand48 == 0)
    {
        seed = 1202107158;
        //seed = (unsigned)time(NULL);
        srand48_r (seed, &drand_buf);
        iset_drand48 = 1; printf("iset_drand48=%d\n",iset_drand48);
    }

    for(i=0;i<nrand;i++)
        drand48_r (&drand_buf, &random[i]);
    return random;
}

int **drand48_int2d(int n1,int n2,int range)
{
    int i,j,seed,**random; double tmp;
    static struct drand48_data drand_buf;
    static int iset_drand48=0;

    random = I_malloc2(n1,n2);
    if(iset_drand48 == 0)
    {
        //seed = 1202107158;
        seed = (unsigned)time(NULL);
        srand48_r (seed, &drand_buf);
        iset_drand48 = 1; printf("iset_drand48=%d\n",iset_drand48);
    }

    for(i=0;i<n1;i++)
        for(j=0;j<n2;j++)
        {
            drand48_r (&drand_buf, &tmp);
            random[i][j] = floor(tmp*(range-1));
        }
    return random;
}

double *return_doublerand(int n,double xmin,double xmax)
{
    int i;double *rand;
    rand = return_drand48(n);
    for(i=0;i<n;i++)
        rand[i] = xmin+(xmax-xmin)*rand[i];
    return rand;
}

//double **return_drand48_2d(int n1,int n2)
//{
//    int i,j,seed;
//    static struct drand48_data drand_buf;
//    double **random;
//    static int iset_drand48=0;
//
//    random = D_malloc2(n1,n2);
//    if(iset_drand48 == 0)
//    {
//        //seed = 1202107158;
//        seed = (unsigned)time(NULL);
//        srand48_r (seed, &drand_buf);
//        iset_drand48 = 1; printf("iset_drand48=%d\n",iset_drand48);
//    }
//
//    for(i=0;i<n1;i++)
//        for(j=0;j<n2;j++)
//            drand48_r (&drand_buf, &random[i][j]);
//    return random;
//}


void rand_point2d(double (*group_pos)[5],int ngroup,double xmin,double xmax,double ymin,double ymax)
{
    int i;double *randx,*randy;
    randx = return_doublerand(ngroup,xmin,xmax);
    randy = return_doublerand(ngroup,ymin,ymax);
    for(i=0;i<ngroup;i++)
    {
        group_pos[i][0] = randx[i];
        group_pos[i][1] = randy[i];
    }
    free(randx);
    free(randy);
}

