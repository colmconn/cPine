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

MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

ctrlSubjects="$( cat $CONFIG_DATA/control.subjectList.txt )"
mddSubjects="$( cat $CONFIG_DATA/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"

function makeAutocorrelatedBuckets {
    local groupName=$1
    local subjectList="$2"

## Autocorrelation based seeds
    csvFilePrefix=subjectOrder.$groupName.REML
    noDataCsvFile=$GROUP_DATA/nodata.$groupName.REML.txt
    rm -f $noDataCsvFile
    touch $noDataCsvFile
    
    echo "subject" > $GROUP_DATA/$csvFilePrefix.csv
    bucketListFile=$GROUP_DATA/bucketList.$groupName.REML.txt
    rm -f $bucketListFile
    touch $bucketListFile
    
    for subject in $subjectList; do
	zScoreFile=$DATA/$subject/functional/rsfc/${seedName}_REML.z-score+tlrc.HEAD
	if [ -f $zScoreFile ] ; then
	    echo "$zScoreFile" >> $bucketListFile
	    echo "$subject" >> $GROUP_DATA/$csvFilePrefix.csv
	else 
	    echo "$subject" >> $noDataCsvFile
	    echo "*** WARNING $zScoreFile does not exist!"
	fi
    done ## end of for subject in $subjects; do
    restingStateBucketFile=$GROUP_DATA/restingstate.bucket.$groupName.REML+tlrc
    echo "Making resting state bucket file $restingStateBucketFile"
    
    rm -f ${restingStateBucketFile}.*
    3dbucket -fbuc -prefix $restingStateBucketFile filelist:$bucketListFile
    rm -f $bucketListFile
    
    rm -f ${restingStateBucketFile%%+*}.masked*
    3dcalc -datum float -a $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz -b ${restingStateBucketFile} -expr "a*b" -prefix ${restingStateBucketFile%%+*}.masked
    rm -f ${restingStateBucketFile}.*
    
}

##makeNonAutocorrelatedBuckets "mddAndCtrl" "$subjects"
makeAutocorrelatedBuckets    "mddAndCtrl" "$subjects"

##makeNonAutocorrelatedBuckets "ctrlOnly" "$ctrlSubjects"
makeAutocorrelatedBuckets    "ctrlOnly" "$ctrlSubjects"

##makeNonAutocorrelatedBuckets "mddOnly"  "$mddSubjects"
makeAutocorrelatedBuckets    "mddOnly"  "$mddSubjects"