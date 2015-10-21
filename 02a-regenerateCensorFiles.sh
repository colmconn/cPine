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

task=pine
censorThreshold=0.2

GETOPT_OPTIONS=$( $GETOPT  -o "s:d:c:" --longoptions "subject:,drop:,censorThreshold:" -n ${programName} -- "$@" )
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
#	-p|--polort)
#	    polort=$2; shift 2 ;;
	# -f|--fwhm)
	#     fwhm=$2 ;
	#     sigma=$( echo "scale=10 ; ${fwhm}/2.3548" | bc ) ;
	#     shift 2 ;;
	-c|--censorThreshold)
	    censorThreshold=$2; shift 2 ;;
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

echo "Using $censorThreshold as threshold of censoring motion"

prefixPrefix=${subjectNumber}.pine

cd  $DATA/$subjectNumber/functional

if [ -f $subjectNumber.pine_tsh_vr_motion.1D ] ; then 
    motionFile=$subjectNumber.pine_tsh_vr_motion.1D
else
    motionFile=$subjectNumber.pine_vr_motion.1D
fi

tr=2000

if [[ $drop != 0 ]] ; then 
    echo "*** Dropping first TRs"
    tail -n +$( expr $drop + 1 ) $motionFile > ${prefixPrefix}.motion.1D
    tail -n +$( expr $drop + 1 ) rm.outliers.censor.1D > rm.outliers.censor.trimmed.1D
else
    echo "*** Not dropping any volumes as drop = $drop"
    cat $motionFile > ${prefixPrefix}.motion.1D
fi

# compute de-meaned motion parameters (for use in regression)
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -demean -write motion_demean.1D

# compute motion parameter derivatives (for use in regression)
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -derivative -demean -write motion_deriv.1D

# create censor file motion_${subj}_censor.1D, for censoring motion 
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -set_tr $tr -show_censor_count -censor_prev_TR                           \
    -censor_motion ${censorThreshold} motion_${subjectNumber}

if [[ $drop != 0 ]] ; then
## now combine both the motion and outlier censor files into one
    1deval -a motion_${subjectNumber}_censor.1D -b rm.outliers.censor.trimmed.1D  \
	-expr "a*b" > censor_${subjectNumber}_combined_2.1D
else 
    1deval -a motion_${subjectNumber}_censor.1D -b rm.outliers.censor.1D  \
	-expr "a*b" > censor_${subjectNumber}_combined_2.1D
fi

motionThresholdPrecentage=0.2
numberOfCensoredVolumes=$( cat censor_${subjectNumber}_combined_2.1D | gawk '{a+=(1-$0)}END{print a}' )
totalNumberOfVolumes=$( cat censor_${subjectNumber}_combined_2.1D | wc -l )
twentyPercent=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

if [[ $numberOfCensoredVolumes -gt $twentyPercent ]] ; then 
    echo "A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than 20% (n=$twentyPercent) of all total volumes of this subject" >  00_DO_NOT_ANALYSE_${subjectNumber}.txt
    echo "*** WARNING: This subject will not be analysed due to having more that 20% of their volumes censored."
fi