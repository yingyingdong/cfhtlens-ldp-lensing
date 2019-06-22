double *zbin_lg(int nzbinlg,double zmin,double zmax,double *zstep);
double gain_zbin_lg(double *DAzlg,double z,double zmin,double zstep);
double *DA_zlg(int nzbinlg,double *zbinlg);
double DA_z(double z);
double D_COM(double z);
double DA_z12(double result1,double result2,double z1);
double sigma(double z1,double z2,double zmin,double zstep,double DA1,double *DAz);
