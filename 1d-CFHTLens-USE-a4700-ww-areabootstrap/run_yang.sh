icc mymath.c linklist.c read_yang.c e_profile_calibrationw2_lg_yang.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#readme : linklist=0 hasn't been deleted
#icc mymath.c linklist.c e_profile_calibrationw2_yang_xy.c -lm -I. -lfftw3 -openmp -o read -g    #bin in dx and dy
#icc mymath.c linklist.c e_profile_calibrationw2_lg_yang_xy.c -lm -I. -lfftw3 -openmp -o read -g    #bin in dx and dy
export OMP_NUM_THREADS=5
