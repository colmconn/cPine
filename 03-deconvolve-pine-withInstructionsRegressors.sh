#!/bin/bash

set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

PINE_REGRESSORS_ROOT=$DATA/regressors
task=pine

GETOPT_OPTIONS=$( $GETOPT  -o "s:p:" --longoptions "subject:,polort:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

drop=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-d|--drop)
	    drop=$2; shift 2 ;;
	-p|--polort)
	    polort=$2; shift 2 ;;
#	-f|--fwhm)
#	    fwhm=$2 ;
#	    sigma=$( echo "scale=10 ; ${fwhm}/2.3548" | bc ) ;
#	    shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided."
    exit
fi

if [ -z $polort ] ; then 
    echo "*** ERROR: The polort was not provided."
    exit
fi

prefixPrefix=${subjectNumber}.pine

cd $DATA/$subjectNumber/functional

##if [ 1 == 0 ] ; then

if [ ! -f ${subjectNumber}.pine.preprocessed.MNI+tlrc.HEAD   ] ; then 
    echo "*** $DATA/$subjectNumber/functional/${subjectNumber}.pine.preprocessed.MNI+tlrc.HEAD does not exist. Cannot continue. Exiting."
fi

if [ ! -f $DATA/regressors/${subjectNumber}.acquisition.regressors.tab ] ; then 
    echo "*** $DATA/regressors/${subjectNumber}.acquisition.regressors.tab (so derivative of it cannot exist) does not exist. Cannot continue. Exiting."
fi


goforit=""
for f in ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.normalized.1D \
    ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.normalized.1D  ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.normalized.1D   ; do

    if [ ! -f $f ] ; then 
	echo "*** Cannot find $f. Exiting"
	exit
    fi
            # this replaces all the newlines in a file (containing one
            # entry per line) with a + and uses bc to sum the numbers
            # the sed is to replace 'e' scientific notation with the *10^ notation preferred by bc
    sum=$( ( tr '\n' +; echo 0 ) < $f | sed 's/[eE]/\*10\^/g' | bc  )
            #echo "$f: $sum"
            # if the sum is zero we had no responses for a particular
            # regressor and therefore need a GOFORIT argument to
            # 3dDeconvolve to prevent it dying. Note that we can stop
            # as soon as we find the first regressor that sums to 0
    if [ $sum == 0 ] ; then
        echo "*** WARNING found a regressor ($f) that sums to zero. Will pass -GOFORIT 99 to 3dDeconvolve"
        goforit="-GOFORIT 99"
        break;
    fi
done

3dDeconvolve -x1D_stop $JOBS \
    $goforit \
    -fout -tout \
    -bout -nofull_first \
    -float \
    -censor censor_${subjectNumber}_combined_2.1D \
    -input  ${prefixPrefix}.preprocessed.MNI+tlrc \
    -bucket ${prefixPrefix}.dec.emotive.withInstructions \
    -errts  ${prefixPrefix}.errts.emotive.withInstructions \
    -fitts  ${prefixPrefix}.fitts.emotive.withInstructions \
    -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
    -polort $polort \
    -xjpeg $subjectNumber.matrix.emotive.withInstructions \
    -x1D $subjectNumber.emotive.withInstructions \
    -num_stimts 14 \
    -stim_file  1 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrA.binary.1D  \
    -stim_label 1 instrA \
    -stim_file  2 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrB.binary.1D  \
    -stim_label 2 instrB \
    -stim_file  3 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrC.binary.1D  \
    -stim_label 3 instrC \
    -stim_file  4 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrD.binary.1D  \
    -stim_label 4 instrD \
    -stim_file  5 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.normalized.1D  \
    -stim_label 5 happy \
    -stim_file  6 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.normalized.1D  \
    -stim_label 6 fearful \
    -stim_file  7 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.normalized.1D  \
    -stim_label 7 neutral \
    -stim_file  8 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.normalized.1D  \
    -stim_label 8 sad \
    -stim_file  9 motion_demean.1D\[0\] \
    -stim_base  9 \
    -stim_label 9 Roll \
    -stim_file  10 motion_demean.1D\[1\] \
    -stim_base  10 \
    -stim_label 10 Pitch \
    -stim_file  11 motion_demean.1D\[2\] \
    -stim_base  11 \
    -stim_label 11 Yaw \
    -stim_file  12 motion_demean.1D\[3\] \
    -stim_base  12 \
    -stim_label 12 dS \
    -stim_file  13 motion_demean.1D\[4\] \
    -stim_base  13 \
    -stim_label 13 dL \
    -stim_file  14 motion_demean.1D\[5\] \
    -stim_base  14 \
    -stim_label 14 dP \
    -num_glt 13 \
    -gltsym 'SYM: +fearful -happy' \
    -glt_label 1 'fearfulVsHappy' \
    -gltsym 'SYM: +fearful -neutral' \
    -glt_label 2 'fearfulVsNeutral' \
    -gltsym 'SYM: +fearful -sad' \
    -glt_label 3 'fearfulVsSad' \
    -gltsym 'SYM: +happy -neutral' \
    -glt_label 4 'happyVsNeutral' \
    -gltsym 'SYM: +happy -sad' \
    -glt_label 5 'happyVsSad' \
    -gltsym 'SYM: +neutral -sad' \
    -glt_label 6 'neutralVsSad' \
    -gltsym 'SYM: +0.33*happy +0.33*fearful +0.33*sad -neutral' \
    -glt_label 7 'allEmotiveVsNeutral' \
    -gltsym 'SYM: +instrA -instrB' \
    -glt_label 8 'instrAVsinstrB' \
    -gltsym 'SYM: +instrA -instrC' \
    -glt_label 9 'instrAVsinstrC' \
    -gltsym 'SYM: +instrA -instrD' \
    -glt_label 10 'instrAVsinstrD' \
    -gltsym 'SYM: +instrB -instrC' \
    -glt_label 11 'instrBVsinstrC' \
    -gltsym 'SYM: +instrB -instrD' \
    -glt_label 12 'instrBVsinstrD' \
    -gltsym 'SYM: +instrC -instrD' \
    -glt_label 13 'instrCVsinstrD' \
    

mv -f $subjectNumber.REML_cmd  $subjectNumber.emotive.withInstructions.REML_cmd 
chmod +x $subjectNumber.emotive.withInstructions.REML_cmd 

if [ -f $subjectNumber.emotive.withInstructions.REML_cmd ] ; then
    ./$subjectNumber.emotive.withInstructions.REML_cmd  $goforit
else
    echo "*** $programName ERROR *** Can't find $subjectNumber.emotive.withInstructions.REML_cmd"
fi

if [ -f $subjectNumber.$task.errts.emotive.withInstructions+tlrc.HEAD ] ; then 
    3dFWHMx -combine -detrend -mask ${prefixPrefix}.preprocessed.mask.MNI+tlrc  \
	$subjectNumber.$task.errts.emotive.withInstructions+tlrc.HEAD > $subjectNumber.$task.errts.emotive.withInstructions.blur.est.1D
fi

if [ -f $subjectNumber.$task.errts.emotive.withInstructions_REML+tlrc.HEAD ] ; then 
    3dFWHMx -combine -detrend -mask ${prefixPrefix}.preprocessed.mask.MNI+tlrc  \
	$subjectNumber.$task.errts.emotive.withInstructions_REML+tlrc.HEAD > $subjectNumber.$task.errts.emotive.withInstructions_REML.blur.est.1D
fi


########################################################################################################################################################################################################

goforit=""
for f in \
    ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.not.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.omission.normalized.1D \
    ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.not.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.omission.normalized.1D \
    ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.not.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.omission.normalized.1D \
    ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.not.remembered.normalized.1D ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.omission.normalized.1D \
    ; do

    if [ ! -f $f ] ; then 
	echo "*** Cannot find $f. Exiting"
	exit
    fi
            # this replaces all the newlines in a file (containing one
            # entry per line) with a + and uses bc to sum the numbers
            # the sed is to replace 'e' scientific notation with the *10^ notation preferred by bc
    sum=$( ( tr '\n' +; echo 0 ) < $f | sed 's/[eE]/\*10\^/g' | bc  )
            #echo "$f: $sum"
            # if the sum is zero we had no responses for a particular
            # regressor and therefore need a GOFORIT argument to
            # 3dDeconvolve to prevent it dying. Note that we can stop
            # as soon as we find the first regressor that sums to 0
    if [ $sum == 0 ] ; then
        echo "*** WARNING found a regressor ($f) that sums to zero. Will pass -GOFORIT 99 to 3dDeconvolve"
        goforit="-GOFORIT 99"
        break;
    fi
done

####################################################################################################
########## Including omissions  ####################################################################
####################################################################################################

# 3dDeconvolve $JOBS \
#     $goforit \
#     -fout -tout \
#     -bout -nofull_first \
#     -float \
#     -censor censor_${subjectNumber}_combined_2.1D \
#     -input  ${prefixPrefix}.preprocessed.MNI+tlrc \
#     -bucket ${prefixPrefix}.dec.memory \
#     -errts  ${prefixPrefix}.errts.memory \
#     -fitts  ${prefixPrefix}.fitts.memory \
#     -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
#     -polort $polort \
#     -xjpeg $subjectNumber.matrix.memory \
#     -x1D $subjectNumber.memory \
#     -num_stimts 17 \
#     -stim_file  1 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.remembered.normalized.1D  \
#     -stim_label 1 happyRem \
#     -stim_file  2 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.not.remembered.normalized.1D  \
#     -stim_label 2 happyNotrem \
#     -stim_file  3 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.omission.normalized.1D  \
#     -stim_label 3 happyOmission \
#     -stim_file  4 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.remembered.normalized.1D  \
#     -stim_label 4 fearfulRem \
#     -stim_file  5 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.not.remembered.normalized.1D  \
#     -stim_label 5 fearfulNotrem \
#     -stim_file  6 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.omission.normalized.1D  \
#     -stim_label 6 fearfulOmission \
#     -stim_file  7 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.remembered.normalized.1D  \
#     -stim_label 7 neutralRem \
#     -stim_file  8 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.not.remembered.normalized.1D  \
#     -stim_label 8 neutralNotrem \
#     -stim_file  9 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.omission.normalized.1D  \
#     -stim_label 9 neutralOmission \
#     -stim_file  10 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.remembered.normalized.1D  \
#     -stim_label 10 sadRem \
#     -stim_file  11 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.not.remembered.normalized.1D  \
#     -stim_label 11 sadNotrem \
#     -stim_file  12 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.omission.normalized.1D  \
#     -stim_label 12 sadOmission \
#     -stim_file  13 motion_demean.1D\[1\] \
#     -stim_base  13 \
#     -stim_label 13 Pitch \
#     -stim_file  14 motion_demean.1D\[2\] \
#     -stim_base  14 \
#     -stim_label 14 Yaw \
#     -stim_file  15 motion_demean.1D\[3\] \
#     -stim_base  15 \
#     -stim_label 15 dS \
#     -stim_file  16 motion_demean.1D\[4\] \
#     -stim_base  16 \
#     -stim_label 16 dL \
#     -stim_file  17 motion_demean.1D\[5\] \
#     -stim_base  17 \
#     -stim_label 17 dP \
#     -gltsym 'SYM: +happyRem -happyNotrem' \
#     -glt_label 1 'happyRemVsHappyNotrem' \
#     -gltsym 'SYM: +fearfulRem -fearfulNotrem' \
#     -glt_label 2 'fearfulRemVsFearfulNotrem' \
#     -gltsym 'SYM: +neutralRem -neutralNotrem' \
#     -glt_label 3 'neutralRemVsNeutralNotrem' \
#     -gltsym 'SYM: +sadRem -sadNotrem' \
#     -glt_label 4 'sadRemVsSadNotrem' \
#     -gltsym 'SYM: +0.25*happyRem +0.25*fearfulRem +0.25*neutralRem  +0.25*sadRem -0.25*happyNotrem -0.25*fearfulNotrem -0.25*neutralNotrem  -0.25*sadNotrem' \
#     -glt_label 5 'allRemVsAllNotrem'

3dDeconvolve -x1D_stop $JOBS \
    $goforit \
    -fout -tout \
    -bout -nofull_first \
    -float \
    -censor censor_${subjectNumber}_combined_2.1D \
    -input  ${prefixPrefix}.preprocessed.MNI+tlrc \
    -bucket ${prefixPrefix}.dec.memory.withInstructions \
    -errts  ${prefixPrefix}.errts.memory.withInstructions \
    -fitts  ${prefixPrefix}.fitts.memory.withInstructions \
    -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
    -polort $polort \
    -xjpeg $subjectNumber.matrix.memory.withInstructions \
    -x1D $subjectNumber.memory.withInstructions \
    -num_stimts 18 \
    -stim_file  1 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrA.binary.1D  \
    -stim_label 1 instrA \
    -stim_file  2 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrB.binary.1D  \
    -stim_label 2 instrB \
    -stim_file  3 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrC.binary.1D  \
    -stim_label 3 instrC \
    -stim_file  4 ${PINE_REGRESSORS_ROOT}/$subjectNumber.acquisition.instrD.binary.1D  \
    -stim_label 4 instrD \
    -stim_file  5 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.remembered.normalized.1D  \
    -stim_label 5 happyRem \
    -stim_file  6 ${PINE_REGRESSORS_ROOT}/$subjectNumber.happy.not.remembered.normalized.1D  \
    -stim_label 6 happyNotrem \
    -stim_file  7 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.remembered.normalized.1D  \
    -stim_label 7 fearfulRem \
    -stim_file  8 ${PINE_REGRESSORS_ROOT}/$subjectNumber.fearful.not.remembered.normalized.1D  \
    -stim_label 8 fearfulNotrem \
    -stim_file  9 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.remembered.normalized.1D  \
    -stim_label 9 neutralRem \
    -stim_file  10 ${PINE_REGRESSORS_ROOT}/$subjectNumber.neutral.not.remembered.normalized.1D  \
    -stim_label 10 neutralNotrem \
    -stim_file  11 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.remembered.normalized.1D  \
    -stim_label 11 sadRem \
    -stim_file  12 ${PINE_REGRESSORS_ROOT}/$subjectNumber.sad.not.remembered.normalized.1D  \
    -stim_label 12 sadNotrem \
    -stim_file  13 motion_demean.1D\[0\] \
    -stim_base  13 \
    -stim_label 13 Roll \
    -stim_file  14 motion_demean.1D\[1\] \
    -stim_base  14 \
    -stim_label 14 Pitch \
    -stim_file  15 motion_demean.1D\[2\] \
    -stim_base  15 \
    -stim_label 15 Yaw \
    -stim_file  16 motion_demean.1D\[3\] \
    -stim_base  16 \
    -stim_label 16 dS \
    -stim_file  17 motion_demean.1D\[4\] \
    -stim_base  17 \
    -stim_label 17 dL \
    -stim_file  18 motion_demean.1D\[5\] \
    -stim_base  18 \
    -stim_label 18 dP \
    -num_glt 11 \
    -gltsym 'SYM: +happyRem -happyNotrem' \
    -glt_label 1 'happyRemVsHappyNotrem' \
    -gltsym 'SYM: +fearfulRem -fearfulNotrem' \
    -glt_label 2 'fearfulRemVsFearfulNotrem' \
    -gltsym 'SYM: +neutralRem -neutralNotrem' \
    -glt_label 3 'neutralRemVsNeutralNotrem' \
    -gltsym 'SYM: +sadRem -sadNotrem' \
    -glt_label 4 'sadRemVsSadNotrem' \
    -gltsym 'SYM: +0.25*happyRem +0.25*fearfulRem +0.25*neutralRem  +0.25*sadRem -0.25*happyNotrem -0.25*fearfulNotrem -0.25*neutralNotrem  -0.25*sadNotrem' \
    -glt_label 5 'allRemVsAllNotrem' \
    -gltsym 'SYM: +instrA -instrB' \
    -glt_label 6 'instrAVsinstrB' \
    -gltsym 'SYM: +instrA -instrC' \
    -glt_label 7 'instrAVsinstrC' \
    -gltsym 'SYM: +instrA -instrD' \
    -glt_label 8 'instrAVsinstrD' \
    -gltsym 'SYM: +instrB -instrC' \
    -glt_label 9 'instrBVsinstrC' \
    -gltsym 'SYM: +instrB -instrD' \
    -glt_label 10 'instrBVsinstrD' \
    -gltsym 'SYM: +instrC -instrD' \
    -glt_label 11 'instrCVsinstrD' \


mv -f $subjectNumber.REML_cmd  $subjectNumber.memory.withInstructions.REML_cmd 
chmod +x $subjectNumber.memory.withInstructions.REML_cmd 

if [ -f $subjectNumber.memory.withInstructions.REML_cmd ] ; then
    ./$subjectNumber.memory.withInstructions.REML_cmd  $goforit
else
    echo "*** $programName ERROR *** Can't find $subjectNumber.memory.withInstructions.REML_cmd"
fi

if [ -f $subjectNumber.$task.errts.memory.withInstructions+tlrc.HEAD ] ; then 
    3dFWHMx -combine -detrend -mask ${prefixPrefix}.preprocessed.mask.MNI+tlrc  \
	$subjectNumber.$task.errts.memory.withInstructions+tlrc.HEAD > $subjectNumber.$task.errts.memory.withInstructions.blur.est.1D
fi

if [ -f $subjectNumber.$task.errts.memory.withInstructions_REML+tlrc.HEAD ] ; then 
    3dFWHMx -combine -detrend -mask ${prefixPrefix}.preprocessed.mask.MNI+tlrc  \
	$subjectNumber.$task.errts.memory.withInstructions_REML+tlrc.HEAD > $subjectNumber.$task.errts.memory.withInstructions_REML.blur.est.1D
fi

