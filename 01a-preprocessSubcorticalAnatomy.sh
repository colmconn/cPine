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

GETOPT_OPTIONS=$( $GETOPT  -o "s:,f:" --longoptions "subject::" -n ${programName} -- "$@" )
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
    echo "*** ERROR: The subject ID was not provided."
    exit
fi

if [ -d $DATA/$subjectNumber/anat ] ; then 
    cd $DATA/$subjectNumber/anat
    
    if [ -f $subjectNumber.anat.nii.gz ] ; then 
	[[ ! -d subcortical ]] && mkdir subcortical

	cd subcortical

	[[ ! -f $subjectNumber.anat_struc.std.nii.gz  ]] && ln -sf ../$subjectNumber.anat_struc.std.nii.gz 

	[[ ! -f $subjectNumber.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD    ]] && ln -sf ../$subjectNumber.anat_struc_brain.std.2.MNI.nonlinear+tlrc.HEAD
	[[ ! -f $subjectNumber.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz ]] && ln -sf ../$subjectNumber.anat_struc_brain.std.2.MNI.nonlinear+tlrc.BRIK.gz

	[[ ! -f ${subjectNumber}.subcortical_all_fast_firstseg.nii.gz ]] && run_first_all -v -i $subjectNumber.anat_struc.std.nii.gz -o $subjectNumber.subcortical -a ../$subjectNumber.anat_struc_brain.std.2.MNI.affine.mat

	echo "*** Applying nonlinear warp to segmented data"
	[[ ! -f ${subjectNumber}.subcortical_all_fast_firstseg.2.MNI.nonlinear.nii.gz ]] && applywarp \
	    --ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
	    --in=${subjectNumber}.subcortical_all_fast_firstseg \
	    --warp=../${subjectNumber}.std.2.MNI.warpcoef \
	    --out=${subjectNumber}.subcortical_all_fast_firstseg.2.MNI.nonlinear \
	    --interp=nn
    else 
	echo "*** ERROR: The  $DATA/$subjectNumber/anat/$subjectNumber.anat.nii.gz does not exist. Exiting."
    fi

else
    echo "*** ERROR $DATA/$subjectNumber/anat does not exist. Exiting."
fi
