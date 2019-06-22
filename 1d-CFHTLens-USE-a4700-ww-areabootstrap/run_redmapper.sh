mdir=/home/dfy/local/Romberg
#file=e_profile_calibrationw2_void_jzfield_ad_cra_zbin_dsigma_pz_onetime_omp_r-and-xy #void use
#file=e_profile_singleflag_r-and-xy
file=e_profile_groupflag_r-and-xy
gcc mymath.c linklist.c read_clusters.c get_z_dal.c fft_m.c mask_region.c thread_assign.c output.c ${file}.c $mdir/Romberg.c $mdir/power.c -lm -I. -I $mdir -lgsl -lgslcblas -lfftw3 -fopenmp  -o read -g
