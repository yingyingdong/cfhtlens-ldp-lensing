mdir=/home/dfy/local/Romberg
file=e_profile_cfhtlens_find_voidUSE_star
gcc mymath.c linklist.c read_clusters.c get_z_dal.c fft_m.c mask_region.c thread_assign.c output.c ${file}.c $mdir/Romberg.c $mdir/power.c -lm -I. -I $mdir -lfftw3 -fopenmp  -o ${file} -g
