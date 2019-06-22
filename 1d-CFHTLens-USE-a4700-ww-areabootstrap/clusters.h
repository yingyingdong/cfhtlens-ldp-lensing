double angu_dia_d_ww(double z);
int clusters(double (*pos)[5],int num);
int redmapper(double (*pos)[5],int num);
int read_cluster(int flagc,double (*group_pos)[5]);
SHEARR **shearr_calloc(int n1,int n2);
SHEAR **shear_calloc(int n1,int n2);
double **get_pz(char *name,int nline);
//read_random.c:
int get_randomc_mask(char *code_dir);
int get_randomc_nomask(char *code_dir);
void read_randfore_mask(int ngroup_rand,double (*group_posr)[5],
                char *code_dir);
void read_randfore_nomask(int ngroup_rand,double (*group_posr)[5],
                char *code_dir);
int my_void(double (*pos)[5],int num,int zbin_use,int nxy,double zstep,int nzbin);

