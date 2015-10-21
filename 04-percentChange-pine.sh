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

REML="_REML"

GETOPT_OPTIONS=$( $GETOPT  -o "s:" --longoptions "subject:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

#drop=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	# -d|--drop)
	#     drop=$2; shift 2 ;;
	# -p|--polort)
	#     polort=$2; shift 2 ;;
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

cd $DATA/$subjectNumber/functional

function extractContrast {
    local analysisName="$1"
    local contrasts="$2"

    if [ -f $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD ] ; then 
	for contrast in $contrasts ; do
	    subbrik=$( 3dinfo -label2index "${contrast}#0_Coef" $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD 2> /dev/null )
	    
	    3dcalc -a $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD\[$subbrik\] \
		-expr "a" \
		-prefix $subjectNumber.$task.$contrast.contrast -datum float
	    
	    3drefit -sublabel 0 "${subjectNumber}.${contrast}.contrast" $subjectNumber.$task.$contrast.contrast+tlrc.HEAD
	    
	done
    else
	echo "*** $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD does not exist. Skipping"
    fi
}

function extractPercentChange {
    local analysisName="$1"
    local stimuli="$2"

    if [ -f $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD ] ; then 
	for stimulus in $stimuli ; do
	    subbrik=$( 3dinfo -label2index "${stimulus}#0_Coef" $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD 2> /dev/null )
	    
	    ## We divide by 100 here because the timeseries are grand
	    ## mean scaled to 10,000 which means that the results
	    ## coming out of the deconvolution are 100 times too
	    ## big. This fixes the issue.
	    3dcalc -a $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD\[$subbrik\] \
		-expr "a/100" \
		-prefix $subjectNumber.$task.$stimulus.%cs -datum float
	    
	    3drefit -sublabel 0 "${subjectNumber}.${stimulus}.%cs" $subjectNumber.$task.$stimulus.%cs+tlrc.HEAD
	    
	done
    else
	echo "*** $subjectNumber.$task.dec.$analysisName$REML+tlrc.HEAD does not exist. Skipping"
    fi
}

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad allEmotiveVsNeutral" 
extractContrast "emotive" "$contrasts"

extractPercentChange "emotive" "happy fearful neutral sad" 


contrasts="happyRemVsHappyNotrem fearfulRemVsFearfulNotrem neutralRemVsNeutralNotrem sadRemVsSadNotrem allRemVsAllNotrem"
extractContrast "memory" "$contrasts"

extractPercentChange "memory" "happyRem happyNotrem fearfulRem fearfulNotrem neutralRem neutralNotrem sadRem sadNotrem"
