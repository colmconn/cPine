#!/bin/bash

#set -x 

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

GETOPT_OPTIONS=$( $GETOPT  -o "s:" --longoptions "subject:" -n ${programName} -- "$@" )
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
#	-d|--drop)
#	    drop=$2; shift 2 ;;
#	-p|--polort)
#	    polort=$2; shift 2 ;;
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


prefixPrefix=${subjectNumber}.pine

cd $DATA/$subjectNumber/functional

stimuli="fearful happy neutral sad"
for stimulus in $stimuli ; do
    if [ ! -f $subjectNumber.acquisition.$stimulus.iresp+tlrc.HEAD   ] ; then 
	echo "*** $DATA/$subjectNumber/functional/$subjectNumber.acquisition.happy.iresp+tlrc.HEAD does not exist. Cannot continue. Exiting."
	exit
    fi

    echo "*** Doing NLFIT on $stimulus for subject $subjectNumber"

    3dNLfim \
	-input $subjectNumber.acquisition.$stimulus.iresp+tlrc.HEAD \
	-mask  ${prefixPrefix}.preprocessed.mask.MNI+tlrc \
	-ignore 0 \
	-noise Constant \
	-signal GammaVar \
	-nconstr 0  -1000.0   1000.0 \
	-sconstr 0    0    2 \
	-sconstr 1  -1000.0   1000.0 \
	-sconstr 2     8    9 \
	-sconstr 3     0.15    0.45   \
	-nrand 500 \
	-nbest 10 \
	-rmsmin 1.0 \
	-bucket 0 ${subjectNumber}.gamma.$stimulus \
	-jobs ${CPUS}
done
