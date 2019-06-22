//redmapper.c
//int redmapper(double (*pos)[4],int num,double xrange[2][2]);
//read_yang.c
int yang_cen(double (*pos)[3],int num,double xrange[2][2]);
int yang_group(double (*pos)[3],int num,double xrange[2][2]);
//read_jz.c
LENS *load_jz(char *file ,int nobj,int num,double *zmax,
		        double xmax[2],double ymax[2]);
LENS *load_jz_v1_3run1(char *file ,int nobj,int num,double *zmax,
		        double xmax[2],double ymax[2],int flag,int flagw1);
int w1_galaxy(double (*pos)[3],int ngroup,double xrange[2][2]);
