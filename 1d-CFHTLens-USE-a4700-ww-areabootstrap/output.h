

void mean_shear_r(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
        double gstep[2],int ntr[2],double zbin[2],double rmin,
        double *DAzlg,double zmin,double zsteplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density);

void mean_shear_xy(double cen[4],LENS *Lens,int Nfound,int *found,SHEAR **shear,
        double gstep[2],int ntr[2],double zbin[2],double rr,
        double *DAzlg,double zmin,double steplg,int flagd,int flag_lg,int flag_wt,double pz,int flag_com,int flag_density);

void gain_nr_gstep_o(int ntr[2],double gstep[2],int flag_lg,int flagd,int flag_xy,double *rmin,double *rmax);

void shear_output_r(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],
        double rmin,int flag_lg,int nuse);

void shear_output_xy(int ntr[2],SHEAR **shear,SHEARR *shearr,char *file,double gstep[2],double rr);

void gain_file_name_xy(int flag_xy,int flag_lg,int flagd,int flagc,int wt_sigma,
        int flag_com,char *tmp1,char *tmp2,char *tmpwt,int ngrp_use,
        double mag_limit,double zmin_v,double rcut,int nzbin_v,double zstep);
void gain_file_name_fore(int flag_xy,int flag_lg,int flagd,int flagc,int wt_sigma,
        int flag_com,char *tmp1,char *tmp2,char *tmpwt,int ngrp_use,
        double zmean_v,int nzbin_v,double zstep,double mag1,double mag2,int flag_single);

int my_void_here_o(double (*pos)[6],double rcut,double zstep,double zmin_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);

int my_voidstar_here_o(double (*pos)[6],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);

int my_fore_galaxy_bino(double (*pos)[5],double zstep,double zmean_v,double zmin,double zsteplg,double *DAzlg,double magl1,double magl2);
int my_fore_galaxy_sortbino(double (*pos)[5],double zstep,double zmean_v,double zmin,double zsteplg,double *DAzlg,double magl1,double magl2);

int my_void_meanfind_o(double (*pos)[5],double rcut,double zstep,double zmin_v,int nzbin_v,double mag_limit,
        double mask_ratio,double zmin,double zsteplg,double *DAzlg,double xrange[2][2],double cutratio);

int edge_check(double x,double ex1,double ex2,double y,double ey1,double ey2,double limit);

void good_data_o(double **arr,int nobj,int nfit,LENS *Lens,double *z,double xma[2],double yma[2],double zcut);

LENS *load_data_o(char *file ,int nobj,int num,double **arr,int *nfit,double *zmax,
        double xmax[2],double ymax[2],double zcut,int flag_read);

void save_found_zbin_o(int **Nfound,int nzbin,int nx,int ny,char *name);

void save_void_pos_o(int **Nfound,int **flag_mask,int NDIV[2],int nzbin,double (*group_pos)[6],
        char *name,double zstep,double zmin,double xc[2]);


void linklist_search_sub_pos2_pz_thread(LENS *Lens,double cen[4],SHEAR **shear,double rmax,double rmin,int *hoc,int *ll, \
                double range[2][2],double step[2],double gstep[2],double zbin[2],int ntr[2],double zsteplg,int NDIV[2],int nfit,double *DAzlg,double zmin,int flagd,int flag_lg,int flag_wt,int flag_com,int flag_density,double pz);

void linklist_search_sub_pos2_pz_thread_xy(LENS *Lens,double cen[4],SHEAR **shear,double rmax,double rmin,int *hoc,int *ll, \
                double range[2][2],double step[2],double gstep[2],double zbin[2],int ntr[2],double zsteplg,int NDIV[2],int nfit,double *DAzlg,
                double zmin,int flagd,int flag_lg,int flag_wt,int flag_com,int flag_density,double pz);
FLAG_GROUP *flag_group_malloc_areaboot(int nnboot,int ngroup,double (*pos)[6]);
int my_void_4field_groupflag(double (*pos)[6],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);
int my_random_4field_groupflag(double (*pos)[6],double rcut,double zstep,double zmean_v,int nzbin_v,
        double mag_limit,double mask_ratio,double zmin,double zsteplg,double *DAzlg);
FLAG_GROUP *flag_group_malloc_areajack(int nnboot,int ngroup,int nfield,double (*pos)[6]);

