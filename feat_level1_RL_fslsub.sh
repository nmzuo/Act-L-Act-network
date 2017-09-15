#!/bin/bash

nSubj=30  ## total RL: 467
for sec in `seq 1 1 $nSubj` 
do
    sec1=${sec}
    sec2=`echo "$sec+45-1"|bc -l `
    if ([ $sec2 -gt $nSubj ])
    then
        sec2=$nSubj
    fi
    sec2=${sec}
    echo "nStart:nEnd " $sec1 $sec2
    /mnt/software/fsl5.0/bin/fsl_sub  bash  feat_level1_RL.sh  $sec1 $sec2  
done
echo "Finished"
