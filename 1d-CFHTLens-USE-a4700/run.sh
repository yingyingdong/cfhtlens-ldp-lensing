#icc mymath.c linklist.c e_profile_calibrationw2_lg.c -lm -I. -lfftw3 -openmp -o read -g
#export OMP_NUM_THREADS=32
date
for mag in -21. -21.5 -22
do
    for nfield in 0 1 2 3 
    #for nfield in 2 3 
    do
        ./e_profile_cfhtlens_find_voidUSE_star ${nfield} ${mag} & 
    done
done
date

#date
#for mag in -20.5 -21 -21.5
#do
#    ./read ${mag}
#done
#date

