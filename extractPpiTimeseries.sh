#!/bin/bash

# set -x 

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
CONFIG_DATA=$DATA/config
PPI_SEEDS_DATA=$CONFIG_DATA/ppi_seeds
scriptsDir=${ROOT}/scripts
GroupResultsDir=$ROOT/data/Group.results
PINE_REGRESSORS_ROOT=$DATA/regressors
task=pine

groups="mddAndCtrl"

usedFwhm=4.2

GETOPT_OPTIONS=$( $GETOPT  -o "s:l:p:c:t:" --longoptions "subject:,seedlist:,polort:,contrast:,stimuli:" -n ${programName} -- "$@" )
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
	-p|--polort)
	    polort=$2; shift 2 ;;
	-c|--contrast)
	    contrast=$2; shift 2 ;;
	-t|--stimuli)
	    stimuli=$2; shift ;;
	-l|--seedlist)
	    seedList=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done


if [ -z $stimuli ] ; then 
    echo "No stimuli provided. Trying to figure these out from the contrast ($contrast)."
    stimuli=$( echo $contrast | sed "s/Vs/ /" |  tr "[:upper:]" "[:lower:]" )
    echo "Setting stimuli to be: $stimuli"
fi

if [ ! -f $seedList ] ; then
    echo "*** ERROR: The seed list file does not exit. Exiting"
    exit
else 
    seeds=$( eval echo $( cat $seedList | sed '/^#/d' ) )
fi

echo "*** Computing PPI for the following seeds:"
echo $seeds

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided. Exiting"
    exit
fi


## the subject's directory from teh resting state analysis
RSFC_SUBJECT_DIR=${RSFC_DATA}/${subjectNumber}/

prefixPrefix=${subjectNumber}.pine

cd $DATA/$subjectNumber/functional/ppi

if [ ! -f ${prefixPrefix}.cleaned.errts_REML+tlrc.HEAD ] ; then
    echo "No such file : $DATA/$subjectNumber/functional/ppi/${prefixPrefix}.cleaned.errts_REML+tlrc.HEAD"
    exit
fi

fLabel="ROIinteraction.group.F-value"
for seed in $seeds ; do

    seedName=${seed##*/}
    seedName=${seedName%%+*}
    ii=${seedName##roi}
    
    mask=$seed
    echo "*** 1. Extracting the cleaned timeseries from the seed ROI"
    3dROIstats -1Dformat -quiet -mask ${mask} ${prefixPrefix}.cleaned.errts_REML+tlrc.HEAD \
	> ${prefixPrefix}.${seedName}.seed.${contrast}.contrast.cleaned.seedTimeseries.1D

    ## now work out the name of the clorder files for this seed and contrast combination

    ## this code (in the if statement below) that checks the length of
    ## the ii variable is specifically to handle the case of the a
    ## priori seeds where the trimming of roi from the beginning to
    ## get the seed number will not change the seedName variable,
    ## leaving roi very long rather than just with just one or two
    ## characters that would indicate the numeric index of the
    ## seed. Note that no more than 99 (or 100 if counting from 0)
    ## seeds can be accomodated in this code. But I'd imagine that
    ## that should be more than enough.
    if [ ${#ii} -le 2 ] ; then 
	suffix=fwhm${usedFwhm}.$task.$groups.roi${ii}.seed.${contrast}.$fLabel
    else 
	suffix=fwhm${usedFwhm}.$task.$groups.${ii}.seed.${contrast}.$fLabel
    fi
    clusterFile="$GroupResultsDir/ppi/clorder.$suffix+tlrc.HEAD"
    if [ -f $clusterFile ] ; then 
	echo "Found $clusterFile"
    else 
	echo "Couldn't find $clusterFile"
    fi

    echo "*** 2. Extracting the cleaned timeseries from all of the ROIs connected to the seed according the the PPI analysis"
    3dROIstats -1Dformat -quiet -mask $clusterFile ${prefixPrefix}.cleaned.errts_REML+tlrc.HEAD  \
	> ${prefixPrefix}.${seedName}.seed.${contrast}.contrast.cleaned.clusterTimeseries.1D

done
