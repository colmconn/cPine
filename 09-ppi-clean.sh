#!/bin/bash

set -x 

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


GETOPT_OPTIONS=$( $GETOPT  -o "s:" --longoptions "subject:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided. Exiting"
    exit
fi

## the subject's directory from teh resting state analysis
RSFC_SUBJECT_DIR=${RSFC_DATA}/${subjectNumber}/

prefixPrefix=${subjectNumber}.pine

cd $DATA/$subjectNumber/functional/ppi

if [ ! -f ${prefixPrefix}.cleaned+tlrc.HEAD ] ; then
    echo "No such file : $DATA/$subjectNumber/functional/ppi/${prefixPrefix}.cleaned+tlrc.HEAD"
    exit
fi

## -input ${prefixPrefix}.preprocessed.MNI+tlrc \
3dDeconvolve -x1D_stop \
    ${JOBS} \
    -rout -tout -bout -float \
    -input  ${prefixPrefix}.cleaned+tlrc \
    -bucket ${prefixPrefix}.cleaned.dec \
    -errts  ${prefixPrefix}.cleaned.errts \
    -x1D    ${prefixPrefix}.cleaned.xmat.1D \
    -mask   ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
    -censor censor_${subjectNumber}_combined_2.1D \
    -num_stimts 6 \
    -stim_file  1 motion_demean.1D\[0\] \
    -stim_base  1 \
    -stim_label 1 Roll \
    -stim_file  2 motion_demean.1D\[1\] \
    -stim_base  2 \
    -stim_label 2 Pitch \
    -stim_file  3 motion_demean.1D\[2\] \
    -stim_base  3 \
    -stim_label 3 Yaw \
    -stim_file  4 motion_demean.1D\[3\] \
    -stim_base  4 \
    -stim_label 4 dS \
    -stim_file  5 motion_demean.1D\[4\] \
    -stim_base  5 \
    -stim_label 5 dL \
    -stim_file  6 motion_demean.1D\[5\] \
    -stim_base  6 \
    -stim_label 6 dP

mv -f $subjectNumber.REML_cmd $subjectNumber.cleaned.REML_cmd 
chmod +x $subjectNumber.cleaned.REML_cmd
./$subjectNumber.cleaned.REML_cmd
