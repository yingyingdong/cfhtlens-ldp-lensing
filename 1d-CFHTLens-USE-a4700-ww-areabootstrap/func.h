#define pi 3.1415926
//#define nboot 100 //for ggsort
#define nboot 400//100
#define nregion 300485//90000//300485//9207
#define griddec 50//50 
#define gridra 50//50 
#define stepdec 16
#define stepra 16 
void myfree(void *mem);
FILE *logfile;
extern FILE *logfile;
#define myfopen(filepointer,filename,filemode) if(!((filepointer)=fopen(filename,filemode))){ fprintf(logfile,"Error opening file '%s'\n",filename);	fflush(logfile); exit(1);	}

extern void load_source_data(int npic,int nstep,float (**position)[2],char *name,int flag);
void load_source_data1(int npic,int nstep,float ***position,char *name,int flag);

float *matrix_new(float kapa,float g1,float g2,float x[2]);
float *matrix_NEW(float kapa,float gama1,float gama2,float *x);
void matrix_2(float kapa,float g1,float g2,float *x);
void matrix(float kapa,float g1,float g2,float x[2]);
void save_shear_data(int npic,int nstep,float (**position)[2],char *name,int flag);
void save_shear_3d_data(int npic,int nstep,float ***position,char *name,int flag);

//read fits file

float *Readimage(char *filename,int *bufsize,int *npix1,int *npix2,float *min,float *max);
void printerror( int status);

float abs_max_pos(float (*pos)[2],int nstep);
float **conv(float **x,float **y,int N1,int N2,int M1,int M2);
void *conv_2d_1(void *x,float *y,int N1,int N2,int M1,int M2,int flag);
void *conv_2d_1_keep(void *x,void *y,int N1,int N2,int M1,int M2,int flag);

//float **float_malloc2(int n1,int n2);
float ***float_malloc3(int n1,int n2,int n3);
void free_float_malloc3(float ***arr,int n1,int n2);
void gaussian_2d(int nsize,float sigma1,float sigma2,float **gauss,int ncx,int ncy);
void gaussian_2d_1(int nsize,float sigma1,float sigma2,void *gauss,int ncx,int ncy,int flag);
void multisize_2d(float **arr,int nsize,int n_ori);
void *multisize_2d_1(void *arr,int nsize,int n_ori,int flag);
void equal_1darr(void *arr1,void *arr2,int n,int flag);
void equal_2darr(float **arr1,float **arr2,int nx,int ny);

//void real_complex_2d_1(fftw_complex *a1,float *a2,int n);
float sum_2d_real(void **arr,int nx,int ny,int flag);
double sum_1d_real(void *arr,int n,int flag);


//mymath.c
void *mymalloc(size_t n);
void *mycalloc(size_t n,int i);
void *mymalloc_complex(size_t n);
void myfree(void *mem);
void *mymalloc_C(size_t n);
int **I_malloc2(int n1,int n2);
int ***I_malloc3(int n1,int n2,int n3);
int **I_ones2(int n1,int n2);
int *I_ones(int n1);
int **I_calloc2(int n1,int n2);
int ***I_calloc3(int n1,int n2,int n3);
float **F_malloc2(int n1,int n2);
float **F_calloc2(int n1,int n2);
double **D_calloc2(int n1,int n2);
double **D_malloc2(int n1,int n2);
fftw_complex **C_malloc2(int n1,int n2);
void free_I_malloc2(int **arr,int n1);
void free_D_malloc2(double **arr,int n1);
void free_F_malloc2(float **arr,int n1);
double ***D_malloc3(int n1,int n2,int n3);
float distance_2dF(float x[2],float y[2]);
double distance_2dD(double x[2],double y[2]);
//void free_C_malloc2(fftw_complex **arr,int n1);
int sum_2d_real_I(int **arr,int nx,int ny);
double sum_2d_real_D(double **arr,int nx,int ny);
float sum_2d_real_F(float **arr,int nx,int ny);
double mean_2d_real_D(double **arr,int nx,int ny);
void initial_1d(void *arr,int n,int flag);
int MAX_I(int a,int b);
int MIN_I(int a,int b);
double MIN_D(double a,double b);
float MAX_F(float a,float b);
double MAX_D(double a,double b);
int *sort_shift(int nsize);
void shift_1d_F(float *arr,int nsize,int *sort);
//void shift_1d_C(fftw_complex *arr,int nsize,int *sort);
float rotation_x(float theta,float x,float y);
float rotation_y(float theta,float x,float y);
void save_int2d(int **arr,int nx,int ny,char *name);
int **load_int2d(int nx,int ny,char *name);
void save_double2d(double **arr,int nx,int ny,char *name);
void save_double2d4(double (*arr)[4],int nx,int ny,char *name);
double **load_double2d(int nx,int ny,char *name);
void GET_line(char *file,int *nobj);
double *return_drand48(int nrand);
void rand_point2d(double (*group_pos)[5],int ngroup,double xmin,double xmax,double ymin,double ymax);
int **drand48_int2d(int n1,int n2,int range);
