#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine/}
DATA=$ROOT/data
GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

logDir=${DATA}/log

task=pine

ctrlSubjects="$( cat ../data/config/control.subjectList.txt )"
mddSubjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad allEmotiveVsNeutral"

for contrast in $contrasts ; do

    maskFile=$GROUP_RESULTS/clorder.fwhm4.2.pine.mddAndCtrl.${contrast}.group.F-value+tlrc.HEAD

    if [[ $contrast == "fearfulVsHappy" ]] ; then
	stimuli="fearful happy"
    elif [[ $contrast == "fearfulVsNeutral" ]] ; then
	stimuli="fearful neutral"
    elif [[ $contrast == "fearfulVsSad" ]] ; then
	stimuli="fearful sad"
    elif [[ $contrast == "happyVsNeutral" ]] ; then
	stimuli="happy neutral"
    elif [[ $contrast == "happyVsSad" ]] ; then
	stimuli="happy sad"
    elif [[ $contrast == "neutralVsSad" ]] ; then
	stimuli="neutral sad"
    elif [[ $contrast == "allEmotiveVsNeutral" ]] ; then
	stimuli="happy fearful neutral sad"
    fi

    for subject in $subjects ; do 
	if [ ! -f $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then 
	    for stimulus in $stimuli ; do 
		roistatsFile=roistats.${subject}.acquisition.${stimulus}.stimulus.${contrast}.contrast.iresp.txt
		( cd $DATA/$subject/functional; 3dROIstats  -nobriklab -mask $maskFile ${subject}.acquisition.${stimulus}.iresp+tlrc.HEAD > $roistatsFile )
	    done
	fi
    done
done