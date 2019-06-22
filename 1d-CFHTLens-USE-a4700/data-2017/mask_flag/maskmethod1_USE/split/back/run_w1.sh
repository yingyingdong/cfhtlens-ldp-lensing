zmean=0.45
for ifield in 1
do
    for dz in 0.2
    do
        for mag in -21 
        do
            for rmax in 0.0166667 0.025
            do
                python split-void-w1.py $ifield $zmean $mag $rmax
            done
        done
    done
done
