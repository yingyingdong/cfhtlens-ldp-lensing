int *thread_assign(int nthread,int n);
int *thread_assign_fixnum(int num,int n,int *nthr);
SHEAR ****shear_thread_assign_nboot(int nthread,int ntr[2]);
SHEARR ***shearr_thread_assign_nboot(int nthread,int ntr[2]);
void zero_shear_thread(SHEAR ***shear,SHEARR **shearr,int **nuse,int ntr[2],int nthread);
void zero_shear_ntr(SHEAR **shear,int ntr[2]);
void zero_nfind_thread(int ***nfind,int n1,int ntr[2]);
void zero_shear_thread_nboot(SHEAR ****shear,SHEARR ***shearr,int **nuse,int ntr[2],int nthread);
void add_nfind_thread(int ****nfind,int nthread,int nzbinvoid,int ntr[2]);
void add_nuse_thread(int **nuse,int nthread,int iboot);
void add_shear_thread(int nthread,SHEAR ****shear,int ntr[2],int iboot);
void add_shear_to_iboot(SHEAR **tmp_shear,SHEAR ***shear,FLAG_GROUP *flag_group,int *nuse,int ntr[2],int index);
