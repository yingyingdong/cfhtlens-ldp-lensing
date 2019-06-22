typedef struct
{
	double e1;
	double e2;
	double sn;
	double de;
	double mag;
	double flag;
	double angle;
	double pos[3];
}LENS;


typedef struct
{
	int nrr;
	double gt;
	double gr;
	double de;
	double err1;
	double err2;
	double err3;
	double dx;
	double dy;
	double dr;
	double angle;
	double wtc;
}SHEAR;

typedef struct
{
	int nrr;
	double gt;
	double gr;
	double de;
	double err1;
	double err2;
	double dr;
	double angle;
}SHEARR;




//linklist:
void makell_sub(LENS *Lens,int nfit,int *hoc,int *ll,int NDIV[2],double range[2][2],\
		        double step[2],double **xgrid);
int *linklist_search_sub(LENS *Lens,int cenid,double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
		            double range[2][2],double step[2],int NDIV[2],int nfit);
int *linklist_search_sub_pos(LENS *Lens,double cen[2],double rmax,double rmin,int *Nfind,int *hoc,int *ll, \
		                double range[2][2],double step[2],int NDIV[2],int nfit);
int *linklist_search_sub_pos2(LENS *Lens,double cen[2],double rmax,int *Nfind,int *hoc,int *ll, \
		                double range[2][2],double step[2],int NDIV[2],int nfit);
void check_num(LENS *Lens,int NDIV[2],double **xgrid,int *hoc,int *ll);
void check_search_sub(int Nfound,int *found,int nfit);

