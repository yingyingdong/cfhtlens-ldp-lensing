#define psf_ratio 0.5
int ELEM(int r,int c,int N);
fftw_complex *zero_ima(int nx,int ny);
double pow_complex(fftw_complex image);
void pow_2darr(fftw_complex *ima_f,int nx,int ny);
float sum_2d(fftw_complex *arr,int nx,int ny);
double mean_2d_complex(fftw_complex *arr,int nx,int ny,int flag);
fftw_complex *make_2dfft(fftw_complex *ima,int nx,int ny,fftw_plan p,int phase);
fftw_complex *make_2difft(fftw_complex *ima,int nx,int ny,fftw_plan p);
void normal_fftifft(fftw_complex *ima,int nx,int ny);
void normal_ima(fftw_complex *ima,int nx,int ny);
fftw_complex *deconvolve(fftw_complex *image_f,fftw_complex *psf_f,fftw_complex *gauss_f,int nx);
void shift_phase_2d(fftw_complex *data, int n1,int n2);
fftw_complex *fftshift_2d(fftw_complex *data, int nx,int ny);
void write_imagepow_txt(fftw_complex *image,int nx,int ny,char *name);
void gain_gauss_fft(fftw_complex *psf_f,fftw_complex *gal_f,fftw_complex *noise_f,
        int nx,float *kx,float *gamma1,float *gamma2,float *gamma3,float h[2]);
void Moffat_profile(fftw_complex *psf,int npsf,int nx,float rd,float rc,float ncx,float ncy);
float Moffat(float rij,float rd,float rc);
fftw_complex *convolve(fftw_complex *image_f,fftw_complex *psf_f,int nx);
fftw_complex *fft_convolve(fftw_complex *image_f,fftw_complex *psf_f,int nx,int ny);
double possion_thresh(fftw_complex *ima,double thresh,int nx,int ny);
void add_gauss_noise(struct drand48_data *drand_buf,fftw_complex *image,
                fftw_complex *nosie,int nx,float sn);
void save_complex2d(fftw_complex *arr,int nx,int ny,char *name);
