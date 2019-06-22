fftw_complex *gain_gauss(double a,double **xgrid,int nx,int ny,double xc,double yc);
int **mask_region(int nx,int ny,int *hoc,int *ll,double step[2],
        double **xgrid,double xc,double yc,int mask_method,int nfield);
int **load_flag_mask(int nfield,int nx,int ny,int load_mask);
