zmean=0.45
#for ifield in 1 2 3 4
#do
#    for dz in 0.2
#    do
#        for mag in -21 -21.5 -22
#        do
#            for rmax in 0.0166667 0.025
#            do
#                python split-void-groupgrid.py $ifield $zmean $mag $rmax 
#            done
#        done
#    done
#done

for mag in -21 -21.5 -22
do
    for rmax in 0.0166667 0.025
    #for rmax in 0.025
    do
        python changeflag-4w-groupgrid.py $zmean $mag $rmax
    done
done
