#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

[[ ! -d $GROUP_DATA ]] && mkdir $GROUP_DATA

task=pine

ctrlSubjects="$( cat ../data/config/control.subjectList.txt )"
mddSubjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"

function makeContrastBuckets {
    local groupName="$1"
    local contrastsName="$2"
    local contrast="$3"
    local subjectList="$4"

    csvFilePrefix=subjectOrder.$groupName.${contrastsName}.REML
    noDataCsvFile=$GROUP_DATA/nodata.$groupName.${contrastsName}.REML.txt
    doNotAnalyzeCsvFile=$GROUP_DATA/doNotAnalyse.$groupName.${contrastsName}.REML.txt
    rm -f $noDataCsvFile
    touch $noDataCsvFile
    rm -f $doNotAnalyzeCsvFile
    touch $doNotAnalyzeCsvFile
    
    echo "subject" > $GROUP_DATA/$csvFilePrefix.csv
    bucketListFile=$GROUP_DATA/bucketList.$groupName.REML.txt
    rm -f $bucketListFile
    touch $bucketListFile
    
    for subject in $subjectList; do
	contrastFile=$DATA/$subject/functional/${subject}.${task}.${contrast}.withInstructions.contrast+tlrc.HEAD
	if [ ! -f $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then
	    if [ -f $contrastFile ] ; then
		echo "$contrastFile" >> $bucketListFile
		echo "$subject" >> $GROUP_DATA/$csvFilePrefix.csv
	    else 
		echo "$subject" >> $noDataCsvFile
		echo "*** WARNING $contrastFile does not exist!"
	    fi
	else 
	    echo "$subject" >> $doNotAnalyzeCsvFile
	    echo "*** WARNING: Found $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt. Skipping ${subject}."
	fi
    done ## end of for subject in $subjects; do
    
    bucketFile=$GROUP_DATA/${task}.bucket.$groupName.${contrastsName}.REML+tlrc
    echo "Making Pine bucket file $bucketFile"
    
    rm -f ${bucketFile}.*
    3dbucket -fbuc -prefix $bucketFile filelist:$bucketListFile
    #rm -f $bucketListFile
    
    rm -f ${bucketFile%%+*}.masked*
    3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${bucketFile} -expr "a*b" -prefix ${bucketFile%%+*}.masked
    rm -f ${bucketFile}.*
}

function makePercentageChangeBuckets {
    local groupName="$1"
    local percentageChangesName="$2"
    local percentageChange="$3"
    local subjectList="$4"

    csvFilePrefix=subjectOrder.$groupName.${percentageChangesName}.REML
    noDataCsvFile=$GROUP_DATA/nodata.$groupName.${percentageChangesName}.REML.txt
    doNotAnalyzeCsvFile=$GROUP_DATA/doNotAnalyse.$groupName.${percentageChangesName}.REML.txt
    rm -f $noDataCsvFile
    touch $noDataCsvFile
    rm -f $doNotAnalyzeCsvFile
    touch $doNotAnalyzeCsvFile
    
    echo "subject" > $GROUP_DATA/$csvFilePrefix.csv
    bucketListFile=$GROUP_DATA/bucketList.$groupName.REML.txt
    rm -f $bucketListFile
    touch $bucketListFile
    
    for subject in $subjectList; do
	percentageChangeFile=$DATA/$subject/functional/${subject}.${task}.${percentageChange}.withInstructions.%cs+tlrc.HEAD
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
    
    bucketFile=$GROUP_DATA/${task}.bucket.$groupName.${percentageChangesName}.%cs.REML+tlrc
    echo "Making Pine bucket file $bucketFile"
    
    rm -f ${bucketFile}.*
    3dbucket -fbuc -prefix $bucketFile filelist:$bucketListFile
    #rm -f $bucketListFile
    
    rm -f ${bucketFile%%+*}.masked*
    3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${bucketFile} -expr "a*b" -prefix ${bucketFile%%+*}.masked
    rm -f ${bucketFile}.*
}

## emotive contrasts

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad allEmotiveVsNeutral happyRemVsHappyNotrem fearfulRemVsFearfulNotrem neutralRemVsNeutralNotrem sadRemVsSadNotrem allRemVsAllNotrem"

for contrast in $contrasts ; do 

    makeContrastBuckets "mddAndCtrl" "$contrast" "$contrast" "$subjects"
    makeContrastBuckets "ctrlOnly"   "$contrast" "$contrast" "$ctrlSubjects"
    makeContrastBuckets "mddOnly"    "$contrast" "$contrast" "$mddSubjects"

done

percentageChangeList="happy fearful sad neutral happyRem happyNotrem fearfulRem fearfulNotrem sadRem sadNotrem neutralRem neutralNotrem"

for pct in $percentageChangeList ; do 

    makePercentageChangeBuckets "mddAndCtrl" "$pct" "$pct" "$subjects"
    makePercentageChangeBuckets "ctrlOnly"   "$pct" "$pct" "$ctrlSubjects"
    makePercentageChangeBuckets "mddOnly"    "$pct" "$pct" "$mddSubjects"

done
