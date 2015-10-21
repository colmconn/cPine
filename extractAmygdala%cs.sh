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

task=pine

ctrlSubjects="$( cat ../data/config/control.subjectList.txt )"
mddSubjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrlSubjects $mddSubjects"

cd $DATA/$subject/

for subject in $subjects ; do
    cd $DATA/$subject/

    if [ ! -f $DATA/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then

	cd anat/subcortical
	
	## 18 is the left amygdala
	## 54 is the right amygdala 
	## see http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIRST/UserGuide#Labels
 
	[[ ! -f ${subject}.amygdala.mask+tlrc.HEAD ]] && \
	    3dcalc -a ${subject}.subcortical_all_fast_firstseg.2.MNI.nonlinear.nii.gz -expr "equals(a, 18) + 2 * equals(a, 54)" -prefix ${subject}.amygdala.mask

	cd ../../functional
	for stimulus in fearful happy neutral sad ; do
	    if [ -f ${subject}.${task}.${stimulus}.%cs+tlrc.HEAD ] ; then
		3dROIstats -nobriklab -mask ../anat/subcortical/${subject}.amygdala.mask+tlrc.HEAD ${subject}.${task}.${stimulus}.%cs+tlrc.HEAD \
		    > roiStats.${subject}.${stimulus}.amygdala.mask.txt
	    fi
	done
    fi
done