#icc mymath.c linklist_jz.c e_profile_jzhang_Yang_lg_back.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c e_profile_jzhang_Yang_lg_e.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c e_profile_jzhang_Yang_lg_e_bgm.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c read_yang.c e_profile_jzhang_Yang_lg_e.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c read_yang.c e_profile_jzhang_Yang_lg_back2.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
icc mymath.c linklist_jz.c read_yang.c read_jz.c e_profile_jzhang_Yang_lg.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c read_yang.c read_jz.c e_profile_jzhang_Yang_lg_v1.3run1.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c read_yang.c read_jz.c e_profile_jzhang_Yang_lg_v1.2.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c e_profile_jzhang_Yang_lg_e_inverse.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c e_profile_jzhang_Yang_xy.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c e_profile_jzhang_redmapper_xy.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle

#icc mymath.c linklist_jz.c read_yang.c read_jz.c e_profile_jzhang_lg.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#icc mymath.c linklist_jz.c read_yang.c read_jz.c e_profile_jzhangv1.3run1_lg.c -lm -I. -lfftw3 -openmp -o read -g      #bin in dr and angle
#readme : linklist=0 hasn't been deleted
export OMP_NUM_THREADS=5
