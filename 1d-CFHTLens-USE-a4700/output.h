LENS *load_data_o(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
        double xmax[2],double ymax[2]);

void mean_shear_r(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
        double gstep[2],int ntr[2],double zbin[2],double rmin,
        double *DAzlg,double zmin,double zsteplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density);

void mean_shear_xy(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
        double gstep[2],int ntr[2],double zbin[2],double rr,double rmin,
        double *DAzlg,double zmin,double steplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density);

void gain_nr_gstep_o(int ntr[2],double gstep[2],int flag_lg,int flagd,int flag_xy,double *rmin,double *rmax);

void shear_output_r(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],
        double rmin,int flag_lg,int nuse);

void shear_output_xy(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],double rr);

void gain_file_name_xy(int flag_xy,int flag_lg,int flagd,int flagc,int wt_sigma,
        int flag_com,char *tmp1,char *tmp2,char *tmpwt,int ngrp_use,
        double mag_limit,double zmin_v,double rcut,int nzbin_v,double zstep);

int my_void_here_o(double (*pos)[5],double rcut,double zstep,double zmin_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);

int my_void_meanfind_o(double (*pos)[5],double rcut,double zstep,double zmin_v,int nzbin_v,double mag_limit,
        double mask_ratio,double zmin,double zsteplg,double *DAzlg,double xrange[2][2],double cutratio);

int edge_check(double x,double ex1,double ex2,double y,double ey1,double ey2,double limit);

void good_data_o(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2]);


LENS *load_data_o(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
        double xmax[2],double ymax[2]);

void save_found_zbin_o(int **Nfound,int nzbin,int nx,int ny,char *name);

void save_void_pos_o(int **Nfound,int **flag_mask,int NDIV[2],int nzbin,double (*group_pos)[4],
        char *name,double zstep,double zmin,double xc[2]);

void save_void_pos_largegrid_o(int **Nfound,int NDIV[2],int nzbin,double (*group_pos)[4],
        char *name,double zstep,double zmean,double xc[2]);
