#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
RSFC_ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/restingstate}
DATA=$ROOT/data
RSFC_DATA=$RSFC_ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts
GroupResultsDir=$ROOT/data/Group.results
PINE_REGRESSORS_ROOT=$DATA/regressors


GETOPT_OPTIONS=$( $GETOPT  -o "s:l:p:c:t:x:b:" --longoptions "subject:,seedlist:,polort:,contrast:,stimuli:,suffix:,seedsubset:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

drop=0
seedsubset=""

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-p|--polort)
	    polort=$2; shift 2 ;;
	-c|--contrast)
	    contrast=$2; shift 2 ;;
	-t|--stimuli)
	    stimuli=$2; shift 2 ;;
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	-x|--suffix)
	    suffix=$2; shift 2 ;;
	-b|--seedsubset)
	    seedsubset=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $polort ] ; then 
    echo "*** ERROR: The polort was not provided."
    exit
fi

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided. Exiting"
    exit
fi

if [ -z "$stimuli" ] ; then 
    echo "No stimuli provided. Trying to figure these out from the contrast ($contrast)."
    stimuli=$( echo $contrast | sed "s/Vs/ /" |  tr "[:upper:]" "[:lower:]" )
    echo "Setting stimuli to be: $stimuli"
fi

if [ ! -f $seedList ] ; then
    echo "*** ERROR: The seed list file does not exit. Exiting"
    exit
else 
    oldSeeds=( $( eval echo $( cat $seedList | sed '/^#/d' ) ) ) 
    seedCount=${#oldSeeds[@]}
    ##echo $seedCount
    ## arrays in bash are zero-based, but the seed files are
    ## one-based. so shift the seed all up by one so thet seed[1]
    ## corresponds to roi1 and so on
    ## declare -a seeds
    for (( ii=${seedCount} ; ii >= 0 ; ii=ii-1 )) 
    do
	##echo "ii: $ii"
	(( newIndex=ii + 1 ))
	##echo "newIndex: $newIndex"
	seeds[$newIndex]="${oldSeeds[$ii]}"
    done
    seeds[0]=""

    # echo "oldSeeds: "
    # for (( ii=0; ii < ${seedCount} + 1 ; ii=ii+1 )) 
    # do
    # 	echo "$ii: ${oldSeeds[$ii]}"
    # done
    
    # ##seeds=${newSeeds[@]}

    # echo "seeds: "
    # for (( ii=0; ii < ${seedCount} + 1 ; ii=ii+1 )) 
    # do
    # 	echo "$ii: ${seeds[$ii]}"
    # done

fi

if [ "x$seedsubset" == "x" ] ; then 
    #echo "*** Computing PPI for the following seeds:"
    #echo ${seeds[*]}
    ##seedCount=${#seeds[*]}
    seedIndices=$( seq 1 $seedCount )
    #echo $seedIndices


else 
    #echo $seedsubset from ${seeds[*]}
    seedIndices=$seedsubset
    #echo $seedIndices
fi

echo "*** Computing PPI for the following seeds:"
for ii in $seedIndices ; do
    echo "$ii: ${seeds[$ii]}"
done

## the subject's directory from the resting state analysis
RSFC_SUBJECT_DIR=${RSFC_DATA}/${subjectNumber}/

prefixPrefix=${subjectNumber}.pine

cd  $DATA/$subjectNumber/functional

# if [ -d ppi ] ; then 
#     mv -f ppi ppi-resultsWithMddSubjectsOnMeds
# fi

mkdir ppi
cd ppi

##if [ 1 == 0 ] ; then 

## link in the motion regressors from the Pine analysis
#ln -sf ../${prefixPrefix}.motion.1D

for f in ../motion_* ../${subject}*motion.1D ; do 
    echo "** Linking in $f"
    [[ ! -f ${f##*/} ]] && ln -sf $f
done


## link in the preprocessed functional data
echo "** Linking in the functional data"
[[ ! -f ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD ]] && ln -sf ../${prefixPrefix}.preprocessed.MNI+tlrc.HEAD
[[ ! -f ${prefixPrefix}.preprocessed.MNI+tlrc.BRIK.gz ]] && ln -sf ../${prefixPrefix}.preprocessed.MNI+tlrc.BRIK.gz

echo "** Linking in the mask"
[[ ! -f ${prefixPrefix}.preprocessed.mask.MNI+tlrc.HEAD ]] && ln -sf ../${prefixPrefix}.preprocessed.mask.MNI+tlrc.HEAD
[[ ! -f ${prefixPrefix}.preprocessed.mask.MNI+tlrc.BRIK.gz ]] && ln -sf ../${prefixPrefix}.preprocessed.mask.MNI+tlrc.BRIK.gz

echo "*** Linking in censor file"
[[ ! -f censor_${subjectNumber}_combined_2.1D ]] && ln -sf ../censor_${subjectNumber}_combined_2.1D

# ## resample the CSF and eroded WM masks to match the functional dataset
# echo "*** Resampling CSF mask to match functional data"
# [[ ! -f ./csf.masked.func.MNI+tlrc.HEAD ]] && 3dresample -master ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD \
#     -inset ../../anat/csf.masked.MNI+tlrc.HEAD \
#     -prefix ./csf.masked.func.MNI
# echo "*** Resampling WM mask to match funtional data"
# [[ ! -f ./wm.eroded.masked.func.MNI+tlrc.HEAD ]] && 3dresample -master ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD \
#     -inset ../../anat/wm.eroded.masked.MNI+tlrc.HEAD \
#     -prefix ./wm.eroded.masked.func.MNI

# 3drefit -space MNI -view tlrc ./csf.masked.func.MNI+tlrc.HEAD
# 3drefit -space MNI -view tlrc ./wm.eroded.masked.func.MNI+tlrc.HEAD

# # extract the average signal from the CSF and WM masks
# echo "*** Extracting average EPI signal from CSF mask"
# [[ ! -f csf.1D ]] && 3dmaskave -q -mask csf.masked.func.MNI+tlrc.HEAD ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD > csf.1D
# echo "*** Extracting average EPI signal from WM mask"
# # 3dmaskave -q -mask wm.eroded.masked.func.MNI+tlrc.HEAD ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD > wm.1D
 
# [[ ! -f ${prefixPrefix}.WM.ts+tlrc.HEAD ]] && 3dLocalstat -nbhd 'SPHERE(30)' -stat mean -mask wm.eroded.masked.func.MNI+tlrc.HEAD -use_nonmask \
#     -prefix ${prefixPrefix}.WM.ts -datum float -quiet ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD

## now detrend the motion, average CSF and average WM timeseries
echo "*** Detrending motion"
[[ ! -f ${prefixPrefix}.motion.detrended.1D  ]] && 3dDetrend -DAFNI_1D_TRANOUT=YES -normalize -prefix ${prefixPrefix}.motion.detrended.1D -polort ${polort} ${prefixPrefix}.motion.1D\'
## echo "*** Detrending CSF"
# [[ ! -f csf.detrended.1D ]] && 3dDetrend -DAFNI_1D_TRANOUT=YES -normalize -prefix csf.detrended.1D -polort ${polort} csf.1D\'
# echo "*** Detrending white matter"
# #3dDetrend -DAFNI_1D_TRANOUT=YES -normalize -prefix wm.detrended.1D -polort ${polort} wm.1D\'
# [[ ! -f ${prefixPrefix}.WM.ts.detrended+tlrc.HEAD ]] && 3dDetrend -normalize -prefix ${prefixPrefix}.WM.ts.detrended -polort ${polort} ${prefixPrefix}.WM.ts+tlrc

# #echo "*** Generating bandpass regressors"
# # create bandpass regressors (instead of using 3dBandpass, say)
# #1dBport -nodata $nvals $tr -band 0.01 0.1 -invert -nozero > bandpass_rall.1D

# [[ ! -f  ${prefixPrefix}_model+tlrc.HEAD ]] && 3dTfitter -RHS ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD \
#     -LHS ${prefixPrefix}.motion.detrended.1D ${prefixPrefix}.WM.ts.detrended+tlrc csf.detrended.1D \
#     -prefix ${prefixPrefix}_model -label roll pitch yaw x y z LocalWM LatVent -fitts ${prefixPrefix}.fitts -quiet

# rm -f ${prefixPrefix}.cleaned+tlrc.*
# [[ ! -f ${prefixPrefix}.cleaned+tlrc.HEAD ]] && 3dcalc -a ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD -b ${prefixPrefix}.fitts+tlrc \
#     -expr "a-b" -prefix ${prefixPrefix}.cleaned

##fi

[[ ! -f ${prefixPrefix}.cleaned+tlrc.HEAD ]] && 3dFourier -prefix ${prefixPrefix}.cleaned \
     -highpass .009 -lowpass .08 \
     -retrend \
     ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD

nvals=$( 3dnvals ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD )
tr=$( 3dinfo -tr ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD )
echo "*** Found $nvals volumes with a TR of ${tr} sec in the preprocessed dataset"

for seedIndex in $seedIndices ; do

    seed=${seeds[$seedIndex]}
    # if [ -z $seed ] ; then 
    # 	## break out of the loop at the first zero length sees. this
    # 	## condition indicates that we've reached the end of the list
    # 	## of seeds.
    # 	break
    # fi
    seedName=${seed##*/}
    seedName=${seedName%%+*}
    
    echo "Running PPI for seed $seedIndex: $seed"
    
    mask=$seed
    echo "*** 1. Extracting the timeseries from the seed ROI"
    3dROIstats -1Dformat -quiet -mask ${mask} ${prefixPrefix}.cleaned+tlrc.HEAD          > ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.1D
    ## 3dROIstats -1Dformat -quiet -mask ${mask} ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD > ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.1D

    echo "*** 2a. Detrending seed timeseries"
    3dDetrend -polort ${polort} -prefix ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.1D \
    	${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.1D\'

    ### the output of this command is the Seed_ts.1D from the AFNI PPI how-to page
    echo "*** 2b. Transposing detrended seed timeseries"
    1dtranspose ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.1D \
    	${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.1D

    echo "*** 3. Generating HRF"
    waver -TR $tr -GAM -inline 1@1 > hrf.1D

    3dTfitter -RHS ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.1D \
    	-FALTUNG hrf.1D ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.convolved 012 -1

    1dtranspose ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.convolved.1D \
    	${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.convolved.transposed.1D 
    
    ##contrastRegressor=${subjectNumber}.ppi.binary.1D
    
    ## use this regressor for the regressors where stimuli are coded
    ## as 1/-1 e.g., fearful = 1 sad = -1

    arr=($(echo $stimuli | tr " " "\n" ))

    #$( echo ${arr[1]}[0:1] | tr '[:lower;]' '[:upper]' )${arr[1]}[1:
    cr="${arr[0]}Vs$(echo ${arr[1]} | gsed 's/.*/\u&/')"

    contrastRegressor=${subjectNumber}.${cr}.contrast.ppi.binary.1D
	
    echo "*** Creating interaction term"
    1deval -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.convolved.transposed.1D \
    	-b $PINE_REGRESSORS_ROOT/${contrastRegressor} -expr "a*b" > \
    	${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.interactionTerm.1D
    
    waver -TR $tr -GAM -peak 1 -numout $nvals -input ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.interactionTerm.1D \
    	> ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.interactionTerm.wavered.1D
    
    	    ## -input ${prefixPrefix}.preprocessed.MNI+tlrc \
    3dDeconvolve -x1D_stop \
    	${JOBS} \
    	-rout -tout -bout -float \
    	-input ${prefixPrefix}.cleaned+tlrc \
    	-bucket ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec \
    	-x1D ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.xmat.1D \
    	-mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
    	-censor censor_${subjectNumber}_combined_2.1D \
    	-num_stimts 13 \
    	-stim_file  1 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.normalized.1D  \
    	-stim_label 1 happy \
    	-stim_file  2 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.normalized.1D  \
    	-stim_label 2 fearful \
    	-stim_file  3 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.normalized.1D  \
    	-stim_label 3 neutral \
    	-stim_file  4 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.normalized.1D  \
    	-stim_label 4 sad \
    	-stim_file  5 motion_demean.1D\[0\] \
    	-stim_base  5 \
    	-stim_label 5 Roll \
    	-stim_file  6 motion_demean.1D\[1\] \
    	-stim_base  6 \
    	-stim_label 6 Pitch \
    	-stim_file  7 motion_demean.1D\[2\] \
    	-stim_base  7 \
    	-stim_label 7 Yaw \
    	-stim_file  8 motion_demean.1D\[3\] \
    	-stim_base  8 \
    	-stim_label 8 dS \
    	-stim_file  9 motion_demean.1D\[4\] \
    	-stim_base  9 \
    	-stim_label 9 dL \
    	-stim_file  10 motion_demean.1D\[5\] \
    	-stim_base  10 \
    	-stim_label 10 dP \
    	-stim_file  11 $PINE_REGRESSORS_ROOT/${contrastRegressor} \
    	-stim_label 11 ${contrast} \
    	-stim_file  12 ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.1D \
    	-stim_label 12 ROIts.${seedName} \
    	-stim_file  13 ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.interactionTerm.wavered.1D \
    	-stim_label 13 ROIinteraction
    
    mv -f $subjectNumber.REML_cmd $subjectNumber.$contrast.REML_cmd 
    chmod +x $subjectNumber.$contrast.REML_cmd
    ./$subjectNumber.$contrast.REML_cmd
    
    # echo "*** Computing z-scores"
    # for label in ROIx${contrast}${suffix} ; do
    # 	subBrikId=$( 3dinfo -label2index "${label}_R^2" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec+tlrc.HEAD 2> /dev/null )
    # 	3dcalc -datum float -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec+tlrc.HEAD\[$subBrikId\] -expr 'log((a+1)/(a-1))/2' \
    # 	    -prefix ${prefixPrefix}.${seedName}.seed.${label}.z-score
    # done

    echo "*** Computing REML z-scores"
    for label in ROIinteraction ; do
    	rqsSubBrikId=$( 3dinfo -label2index "${label}_R^2" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec_REML+tlrc.HEAD 2> /dev/null )
    	betaSubBrikId=$( 3dinfo -label2index "${label}#0_Coef" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec_REML+tlrc.HEAD 2> /dev/null )
	
    	3dcalc -datum float \
    	    -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec_REML+tlrc.HEAD\[$rqsSubBrikId\] \
    	    -b ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec_REML+tlrc.HEAD\[$betaSubBrikId\] \
    	    -expr "ispositive(b)*sqrt(a)-isnegative(b)*sqrt(a)" \
    	    -prefix  ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.REML.corr
	
    	3dcalc -datum float \
    	    -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.REML.corr+tlrc.HEAD -expr 'log((a+1)/(a-1))/2' \
    	    -prefix ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${label}.REML.z-score
    done
done

    























# 3dDeconvolve $JOBS \
#     -polort $polort \
#     -float \
#     -fout -tout \
#     -censor motion_${subjectNumber}_censor.1D \
#     -input  ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD \
#     -bucket ${prefixPrefix}.dec \
#     -errts  ${prefixPrefix}.errts \
#     -fitts  ${prefixPrefix}.fitts \
#     -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc.HEAD \
#     -xjpeg  ${prefixPrefix}.matrix \
#     -x1D    X.xmat.1D \
#     -x1D_uncensored X.nocensor.xmat.1D \
#     -ortvec bandpass_rall.1D bandpass \
#     -x1D_stop  \
#     -num_stimts  8 \
#     -stim_file   1 ${prefixPrefix}.motion.detrended.1D'[0]' -stim_base  1 -stim_label  1 roll_01  \
#     -stim_file   2 ${prefixPrefix}.motion.detrended.1D'[1]' -stim_base  2 -stim_label  2 pitch_01 \
#     -stim_file   3 ${prefixPrefix}.motion.detrended.1D'[2]' -stim_base  3 -stim_label  3 yaw_01   \
#     -stim_file   4 ${prefixPrefix}.motion.detrended.1D'[3]' -stim_base  4 -stim_label  4 dS_01    \
#     -stim_file   5 ${prefixPrefix}.motion.detrended.1D'[4]' -stim_base  5 -stim_label  5 dL_01    \
#     -stim_file   6 ${prefixPrefix}.motion.detrended.1D'[5]' -stim_base  6 -stim_label  6 dP_01    \
#     -stim_file   7 wm.detrended.1D                          -stim_base  7 -stim_label  7 wm \
#     -stim_file   8 csf.detrended.1D                         -stim_base  8 -stim_label  8 csf

# if [ -f $subjectNumber.REML_cmd ] ; then
#     bash $subjectNumber.REML_cmd
# fi




    # for stimulus in $stimuli ; do
	
    # 	contrastRegressor=${subjectNumber}.$stimulus.ppi.binary.1D
	
    # 	## contrastRegressor=${subjectNumber}.${contrast}${suffix}.contrast.ppi.binary.1D
	
    # 	echo "*** Creating interaction term"
    # 	1deval -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.convolved.transposed.1D \
    # 	    -b $PINE_REGRESSORS_ROOT/${contrastRegressor} -expr "a*b" > \
    # 	    ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.${stimulus}.stimulus.interactionTerm.1D
	
    # 	waver -TR $tr -GAM -peak 1 -numout $nvals -input ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.${stimulus}.stimulus.interactionTerm.1D \
    # 	    > ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.${stimulus}.stimulus.interactionTerm.wavered.1D
	
    # 	    ## -input ${prefixPrefix}.preprocessed.MNI+tlrc \
    # 	3dDeconvolve -x1D_stop \
    # 	    ${JOBS} \
    # 	    -rout -tout -bout -float \
    # 	    -input ${prefixPrefix}.cleaned+tlrc \
    # 	    -bucket ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.dec \
    # 	    -x1D ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.xmat.1D \
    # 	    -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
    # 	    -censor censor_${subjectNumber}_combined_2.1D \
    # 	    -num_stimts 13 \
    # 	    -stim_file  1 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.normalized.1D  \
    # 	    -stim_label 1 happy \
    # 	    -stim_file  2 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.normalized.1D  \
    # 	    -stim_label 2 fearful \
    # 	    -stim_file  3 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.normalized.1D  \
    # 	    -stim_label 3 neutral \
    # 	    -stim_file  4 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.normalized.1D  \
    # 	    -stim_label 4 sad \
    # 	    -stim_file  5 motion_demean.1D\[0\] \
    # 	    -stim_base  5 \
    # 	    -stim_label 5 Roll \
    # 	    -stim_file  6 motion_demean.1D\[1\] \
    # 	    -stim_base  6 \
    # 	    -stim_label 6 Pitch \
    # 	    -stim_file  7 motion_demean.1D\[2\] \
    # 	    -stim_base  7 \
    # 	    -stim_label 7 Yaw \
    # 	    -stim_file  8 motion_demean.1D\[3\] \
    # 	    -stim_base  8 \
    # 	    -stim_label 8 dS \
    # 	    -stim_file  9 motion_demean.1D\[4\] \
    # 	    -stim_base  9 \
    # 	    -stim_label 9 dL \
    # 	    -stim_file  10 motion_demean.1D\[5\] \
    # 	    -stim_base  10 \
    # 	    -stim_label 10 dP \
    # 	    -stim_file  11 $PINE_REGRESSORS_ROOT/${contrastRegressor} \
    # 	    -stim_label 11 ${contrast}${suffix} \
    # 	    -stim_file  12 ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.detrended.transposed.1D \
    # 	    -stim_label 12 ROIts.${seedName} \
    # 	    -stim_file  13 ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.contrast.${stimulus}.stimulus.interactionTerm.wavered.1D \
    # 	    -stim_label 13 ROIx${stimulus}
	
    # 	mv -f $subjectNumber.REML_cmd $subjectNumber.$contrast.$stimulus.REML_cmd 
    # 	chmod +x $subjectNumber.$contrast.$stimulus.REML_cmd
    # 	./$subjectNumber.$contrast.$stimulus.REML_cmd
	
    # # echo "*** Computing z-scores"
    # # for label in ROIx${contrast}${suffix} ; do
    # # 	subBrikId=$( 3dinfo -label2index "${label}_R^2" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec+tlrc.HEAD 2> /dev/null )
    # # 	3dcalc -datum float -a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.dec+tlrc.HEAD\[$subBrikId\] -expr 'log((a+1)/(a-1))/2' \
    # # 	    -prefix ${prefixPrefix}.${seedName}.seed.${label}.z-score
    # # done
    # 	echo "*** Computing REML z-scores"
    # 	for label in ROIx${stimulus} ; do
    # 	    rqsSubBrikId=$( 3dinfo -label2index "${label}_R^2" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.dec_REML+tlrc.HEAD 2> /dev/null )
    # 	    betaSubBrikId=$( 3dinfo -label2index "${label}#0_Coef" ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.dec_REML+tlrc.HEAD 2> /dev/null )

    # 	    3dcalc -datum float \
    # 		-a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.dec_REML+tlrc.HEAD\[$rqsSubBrikId\] \
    # 		-b ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.dec_REML+tlrc.HEAD\[$betaSubBrikId\] \
    # 		-expr "ispositive(b)*sqrt(a)-isnegative(b)*sqrt(a)" \
    # 		-prefix  ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.REML.corr

    # 	    3dcalc -datum float \
    # 		-a ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${stimulus}.REML.corr+tlrc.HEAD -expr 'log((a+1)/(a-1))/2' \
    # 		-prefix ${prefixPrefix}.${seedName}.seed.${contrast}${suffix}.${label}.REML.z-score
    # 	done
