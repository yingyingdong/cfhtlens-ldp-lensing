#i=0
for ir in 1 1.5
do
    for zmean in 0.35 0.45 0.512
    do
        for mag in -21 -21.5 -22
        do
            ./read ${mag} ${zmean} ${ir} 
            #((i=$i+1))
        done
    done
done
