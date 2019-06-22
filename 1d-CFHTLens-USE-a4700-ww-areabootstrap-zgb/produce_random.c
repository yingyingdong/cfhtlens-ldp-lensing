

void gain_random_for_void(int **flag_mask,double (*group_pos)[5],int ngroup,double xrange[2][2],
        int NDIV[2],double step[2],char *dir)
{
    int i,j,k,nx[2],mul=1,num;
    double x[2],(*xrand)[4];
    char file[1024];
    xrand = mymalloc(sizeof(double)*4*ngroup*mul);

    srand((unsigned)time(NULL));
    num = 0;
    for(i=0;i<ngroup;i++)
    {
        for(k=0;k<mul;k++)
        {
            do{
                for(j=0;j<2;j++)
                {
                    x[j]  = (rand()/((double)RAND_MAX+1.))*(xrange[j][1]-xrange[j][0])+xrange[j][0];
                    nx[j] = floor((x[j]-xrange[j][0])/step[j]);
                }
            }while(nx[0]>=NDIV[0] && nx[1]>=NDIV[1] && flag_mask[nx[1]][nx[0]]==0);
            for(j=0;j<2;j++) xrand[num][j]=x[j];
            xrand[num][2] = group_pos[i][2];
            num++;
        }
    }

    //sprintf(file,"%s/random_pos_mul%d.dat",dir,mul);
    sprintf(file,"%s/random_pos_mul%d_w1_maskmethod1_nzbin3-step0.006139-zstep0.1-zmin0.2-rmax0.05ANG_magl-20.dat",dir,mul);
    save_double2d4(xrand,3,ngroup*mul,file);
}

