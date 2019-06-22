int *thread_assign(int nthread,int n);
SHEAR ****shear_thread_assign(int nzbin,int nthread,int ntr[2]);
SHEARR ***shearr_thread_assign(int nzbin,int nthread,int ntr[2]);
void zero_shear_thread(SHEAR ****shear,SHEARR ***shearr,int nzbin,int ntr[2],int nthread);
void zero_nfind_thread(int ***nfind,int n1,int ntr[2]);
void add_nfind_thread(int ****nfind,int nthread,int nzbinvoid,int ntr[2]);
void add_nuse_thread(int nuse[],int nthread);
