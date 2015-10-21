#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine/}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
GROUP_RESULTS=$DATA/Group.results
CONFIG_DATA=$DATA/config
PPI_SEEDS_DATA=$CONFIG_DATA/ppi_seeds
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

[[ ! -d $GROUP_DATA ]] && mkdir $GROUP_DATA

task=pine

ctrlSubjects="$( cat ../data/config/control.subjectList.txt )"
mddSubjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"


maskList=""

for subject in $subjects; do
    prefixPrefix=${subject}.${task}
    
    maskFile=$DATA/$subject/functional/${prefixPrefix}.preprocessed.mask.MNI.nii.gz
    if [ ! -f $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then
	maskList="$maskList $maskFile"
    fi
done

3dMean -mask_union -prefix $GROUP_RESULTS/meanMask.union $maskList
 