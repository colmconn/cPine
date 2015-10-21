#!/bin/bash

# set -x 

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

deconvolutionFile=${prefixPrefix}.dec.emotive.tent+tlrc.HEAD

if [ -e $deconvolutionFile ] ; then
    echo "*** Calculating %AUC for subject $subjectNumber"
    
	################
	## define variables for exceptions and the default
	################
    declare -a events
    
    events=( "happy" "fearful" "neutral" "sad" )
    
    for (( ind = 0 ; ind < ${#events[@]} ; ind++ ))
    do
	event=${events[$ind]}
	
	    ## we'll use the variable s for scaling factor
	echo "*** For the $event subbrik"
	gammaFile="${subjectNumber}.gamma.${event}+tlrc.HEAD"
	gammaSubBrikLabel="Signal Area"
	auc=`3dinfo -label2index "$gammaSubBrikLabel" $gammaFile`
	kSubBrikLabel="k"
	k=`3dinfo -label2index "$kSubBrikLabel" $gammaFile`
	
	echo "*** $gammaSubBrikLabel sub-BRIK from GAMMA file is in sub-BRIK: $auc"
	echo "*** $kSubBrikLabel sub-BRIK from GAMMA file is in sub-BRIK: $k"
	
	optev="-k${k} ${gammaFile} -s${auc} ${gammaFile}"
	
	prefix=${prefixPrefix}.${event}.%auc
	
	3dcalc -datum float $optev -expr "(step(k)-step(-k))*100*s/800" -prefix $prefix
	
	3drefit -fim $prefix+tlrc.HEAD
    done
    
else # end of if [ -e $deconvolutionFile ] ; then
    echo "*** SKIPPING calculation of %AUC for subject $subjectNumber: Problems with the timeseries"
    echo "*** $deconvolutionFile does NOT exist"
fi

