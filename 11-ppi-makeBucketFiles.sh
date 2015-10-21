#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
CONFIG_DATA=$DATA/config
PPI_SEEDS_DATA=$CONFIG_DATA/ppi_seeds
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

[[ ! -d $GROUP_DATA ]] && mkdir $GROUP_DATA

task=pine
REML="REML."

GETOPT_OPTIONS=$( $GETOPT  -o "x:" --longoptions "suffix:" -n ${programName} -- "$@" )
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
	-x|--suffix)
	    suffix=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done


ctrlSubjects="$( cat ../data/config/control.subjectList.txt )"
mddSubjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"

#subjects="107_A"

## this file was appropriated from 05-makeBucketFiles.sh hence the
## percentageChange variables sprinkled liberally throughout

function makePercentageChangeBuckets {
    local groupName="$1"
    local percentageChangesName="$2"
    local percentageChange="$3"
    local subjectList="$4"

    csvFilePrefix=subjectOrder.$groupName.${percentageChangesName}.${REML}z-score
    noDataCsvFile=$GROUP_DATA/nodata.$groupName.${percentageChangesName}.${REML}z-score.txt
    doNotAnalyzeCsvFile=$GROUP_DATA/doNotAnalyse.$groupName.${percentageChangesName}.${REML}z-score.txt
    rm -f $noDataCsvFile
    touch $noDataCsvFile
    rm -f $doNotAnalyzeCsvFile
    touch $doNotAnalyzeCsvFile
    
    echo "subject" > $GROUP_DATA/$csvFilePrefix.csv
    bucketListFile=$GROUP_DATA/bucketList.$groupName.${percentageChangesName}.z-score.txt
    rm -f $bucketListFile
    touch $bucketListFile
    
    for subject in $subjectList; do
	percentageChangeFile=$DATA/$subject/functional/ppi/${subject}.${task}.${percentageChange}.${REML}z-score+tlrc.HEAD
	if [ ! -f $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then
	    if [ -f $percentageChangeFile ] ; then
		echo "$percentageChangeFile" >> $bucketListFile
		echo "$subject" >> $GROUP_DATA/$csvFilePrefix.csv
	    else 
		echo "$subject" >> $noDataCsvFile
		echo "*** WARNING $percentageChangeFile does not exist!"
	    fi
	else 
	    echo "$subject" >> $doNotAnalyzeCsvFile
	    echo "*** WARNING: Found $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt. Skipping ${subject}."
	fi
    done ## end of for subject in $subjects; do
    
    bucketFile=$GROUP_DATA/${task}.bucket.$groupName.${percentageChangesName}.${REML}z-score+tlrc
    echo "Making Pine bucket file $bucketFile"
    
    rm -f ${bucketFile}.*
    3dbucket -fbuc -prefix $bucketFile filelist:$bucketListFile
    #rm -f $bucketListFile
    
    rm -f ${bucketFile%%+*}.masked*
    3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${bucketFile} -expr "a*b" -prefix ${bucketFile%%+*}.masked
    ## rm -f ${bucketFile}.*
}

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"

#contrastsRoot="fearfulVsHappy" # fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"

#contrasts="fearfulVsSad"
#contrasts="happyVsNeutral happyVsSad neutralVsSad"

#contrasts="fearfulVsHappyExtractedRightSgAcc" # fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"

#rois="HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala sgacc.left.3mm sgacc.right.3mm"
# anteriorInsula.3mm posteriorInsula.3mm"
#rois="sgacc.left.3mm sgacc.right.3mm"

#rois="apriori.rois.noVmpfc.3mm apriori.rois.withVmpfc.3mm"
rois="apriori.rois.withVmpfc.3mm"

# for cc in $contrastsRoot ; do
#     for rr in $rois ; do 
# 	contrasts="$contrasts $cc$rr"
#     done
# done

## echo $contrasts

for cc in $contrasts ; do
    for rr in $rois ; do
	numberOfSeeds=$( wc -l $PPI_SEEDS_DATA/$cc.$rr.seeds.txt |awk '{print $1}' )

	for (( ii=1; ii <= $numberOfSeeds; ii=ii+1 ))
	do
	    stimuli=$( echo $cc | sed "s/Vs/ /" |  tr "[:upper:]" "[:lower:]" )
	    ##for stimulus in $stimuli ; do
	    for stimulus in ROIinteraction ; do
		##suffix=roi${ii}.seed.${contrast}.ROIx${stimulus}
		suffix=roi${ii}.seed.${cc}${rr}.${stimulus}
		# echo $suffix
		makePercentageChangeBuckets "mddAndCtrl" "$suffix" "$suffix" "$subjects"
	    done
	done
    done
done

exit
## now put together the list of contrasts by a-priori ROIs
if [ "x$rois" == "x" ] ; then 
    echo "*** No a-priori ROIs were defined. Skipping this step."
else 
    for contrast in $contrasts ; do
	for roi in $rois ; do
	    suffix="${roi}.seed.${contrast}.ROIinteraction"
	    makePercentageChangeBuckets "mddAndCtrl" "$suffix" "$suffix" "$subjects"
	done
    done
fi

# percentageChangeList="fearfulVsSad.amygdala.4.seed.ROIxfearfulVsSad"

# for pct in $percentageChangeList ; do 

#     makePercentageChangeBuckets "mddAndCtrl" "$pct" "$pct" "$subjects"
#     makePercentageChangeBuckets "ctrlOnly"   "$pct" "$pct" "$ctrlSubjects"
#     makePercentageChangeBuckets "mddOnly"    "$pct" "$pct" "$mddSubjects"

# done
