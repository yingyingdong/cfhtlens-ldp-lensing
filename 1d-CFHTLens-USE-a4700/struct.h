typedef struct
{
    int gid;
    int flag;
	double e1;
	double e2;
	double M;
	double mag;
	double wt;
	double angle;
	double m;
	double pos[3];
    double pz;
}LENS;


typedef struct
{
	int nrr;
	double gt;
	double gr;
	double wt;
	double err1;
	double err2;
	double dr;
	double angle;
	double m;
}SHEAR;

typedef struct
{
	int nrr;
	double gt;
	double gr;
	double wt;
	double err1;
	double err2;
	double dr;
	double m;
	double angle;
}SHEARR;


#define omegal0 0.732
#define omegam0 0.268
#define H0 72.
#define DH 3000000. //(kpc/h)
#define au (pi/180.)  //error
#define GM 0.00000000430071  //(Mpc/M_sun (km/s)^2)
#define cM 300000. //(km/s)

//linklist:
void makell_sub(LENS *Lens,int nfit,int *hoc,int *ll,int NDIV[2],double range[2][2],\
		double step[2],double **xgrid);
int *linklist_search_sub(LENS *Lens,int cenid,double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
		    double range[2][2],double step[2],int NDIV[2],int nfit);
int *linklist_search_sub_pos(LENS *Lens,double cen[2],double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
		        double range[2][2],double step[2],int NDIV[2],int nfit);
int *linklist_search_sub_pos2(LENS *Lens,double cen[2],double rmax,int *Nfind,int *hoc,int *ll, \
		        double range[2][2],double step[2],int NDIV[2],int nfit);
int *linklist_search_sub_pos2_countzbin12(LENS *Lens,double cen[4],double rmax,int *hoc,int *ll,\
        double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,int **flag_mask,int flagd,double mask_ratio);
int ***linklist_search_sub_pos2_countzbin12_record(LENS *Lens,double cen[5],double rmax,double rmin,double rstep,
        int *hoc,int *ll,double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,
        int **flag_mask,int flagd,int flag_lg,double mask_ratio,int *flag_use,int ntr[2]);
int ***linklist_search_sub_pos2_countzbin12_record_nomask(LENS *Lens,double cen[5],double rmax,double rmin,double rstep,
        int *hoc,int *ll,double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,
        int flagd,int flag_lg,int *flag_use,int ntr[2]);
int *linklist_search_sub_pos2_countzbin12_COM(LENS *Lens,double cen[4],double rmax,int *hoc,int *ll,\
        double range[2][2],double step[2],int NDIV[2],int nzbin,double zstep,double zmean,int **flag_mask,
        double mask_ratio,double DCOM[]);
void check_num(LENS *Lens,int NDIV[2],double **xgrid,int *hoc,int *ll);
void check_search_sub(int Nfound,int *found,int nfit);
