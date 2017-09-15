#!/bin/bash


# generate the confound EVs for CSF and WM


abs () { echo ${1#-};}

# norm_confound in.txt  ## in.txt is the input file
norm_confound() {  
    mymax=0
    mymean=0
    iCount=0
    for i in `cat $1`
    do
        inarr[$iCount]=$i
        mymean=`echo "$mymean  $i" | awk '{printf "%f", $1 + $2}' `
        iCount=$((iCount+1))
    done
    mymean=`echo "$mymean  $iCount" |awk '{printf "%f", $1 / $2}' `

    # get the mean
    for i in $(seq 0 `echo $((iCount-1))`)
    do
        inarr[$i]=`echo "${inarr[$i]} $mymean" |awk '{printf "%f", $1 - $2}'  `
        if (( $(echo " `abs ${inarr[$i]}` > $mymax " |bc -l) )); then
            mymax=`abs ${inarr[$i]}`
        fi
    done

    #output the results
    rm -f $1
    for i in $(seq 0 `echo $((iCount-1))`)
    do
        if (( $(echo " $i < 1 " |bc -l) )); then
            echo "${inarr[$i]} $mymax"| awk '{printf "%f", $1 / $2}' > $1
            echo >> $1
        else
            echo "${inarr[$i]} $mymax"| awk '{printf "%f", $1 / $2}' >> $1
            echo >> $1
        fi
    done
}

featpath='/DATA/239/nmzuo/software/fsl-5.0.9/bin'
csfmask='/DATA/239/nmzuo/software/brat_old/template/fmaskEPI_V2mm_CSF.nii'
wmmask='/DATA/239/nmzuo/software/brat_old/template/fmaskEPI_V2mm_WM.nii'
#tLen=176
orgpath='/DATA/data/HCP/S500/'
newpath='/DATA/239/nmzuo/hcp_S500_defil/'
pubpath='/MNINonLinear/Results/'
#tname=(tfMRI_GAMBLING_RL tfMRI_MOTOR_RL tfMRI_SOCIAL_RL tfMRI_EMOTION_RL tfMRI_LANGUAGE_RL tfMRI_RELATIONAL_RL tfMRI_WM_RL)
tname=(tfMRI_EMOTION_RL)
level1='hp200_s4_level1_RL'

sed_prefsf(){
    #1 for each task
    #2 for each subject
    pubfsf="$newpath$2${pubpath}$1/$level1/level1.fsf"
    curpath=`pwd`

    if [ ! -d "${newpath}$2${pubpath}$1/$level1" ]; then
        mkdir -p ${newpath}$2${pubpath}$1/$level1
    fi
    cd ${newpath}$2${pubpath}$1/$level1
    fslmeants -i $orgpath$2$pubpath$1/$1".nii.gz" -m $csfmask -o csf.txt
    fslmeants -i $orgpath$2$pubpath$1/$1".nii.gz" -m $wmmask -o wm.txt
    norm_confound csf.txt
    norm_confound wm.txt
    paste -d '  '  $orgpath$2$pubpath$1/Movement_Regressors.txt csf.txt wm.txt > mot_csf_wm.txt

    cd $curpath
}


batch_feat(){
    jCount=0
    for ii in `cat subj_30_HCP.txt`  ####VERY strange, the ii var is changed by sed_prefsf() !!!
    do
        jCount=$((jCount+1))
        if ([[ $jCount -lt 1 ]] || [[ $jCount -gt 30 ]])
        then
            continue;
        else
            echo $ii
            for j in ${tname[*]}
            do
    #            sed_prefsf $j $ii # This line could call wrapper bash file, *_fslsub.sh
                curpath=`pwd`
                cd ${newpath}$ii/${pubpath}$j/$level1/
                cat ${newpath}"100307/"${pubpath}$j/$level1/level1.fsf |sed s/100307/$ii/g >  level1.fsf
#                rm -fr level1.feat
                mv level1.feat  level1_biggerArea.feat
                rm -fr ../hp200_s4_level1_RL.feat
                $featpath/feat level1.fsf
                cd $curpath
                echo -e "$j   \c"
            done
            echo ' '
        fi
    done
}

batch_feat 



