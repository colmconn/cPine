#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
DATA=$ROOT/data
GROUP_DATA=$DATA/Group.data
ROIS_DIR=$DATA/rois
GROUP_RESULTS=$DATA/Group.results
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts
logDir=${DATA}/log

task=pine

contrasts="allEmotions"

cd $GROUP_RESULTS
groups="mddAndCtrl"

function pickLatestBucketFile {
    local contrast=$1
    ##restingstate.$groups.$contrast.REML.lme.bucket.*+tlrc.HEAD

    latest=$( ls -1t $task.lme.bucket.%auc.$groups.$contrast.*+tlrc.HEAD | head -1 ) 
    if [ "x$latest" == "x" ] ; then
	exit
    fi
    echo $latest
}

function extractFStatpars {
    local bucket=$1
    local subbrikId=$2

    a=$(3dAttribute BRICK_STATSYM $bucket\[$subbrikId\] )
    b=${a##*(}
    c=${b%%)*}

    echo $( echo $c | tr "," " " )
}

function extractNVoxels {
    local nn=$1
    # corrected p-value
    local cPValue=$2
    ## voxelwise p-value
    local vPValue=$3
    local clustsimPrefix=$4

    #0.100  0.090  0.080  0.070  0.060  0.050  0.040  0.030  0.020  0.010
    if [ $cPValue == 0.100 ] ; then
	pvalueColumn=2
    elif [ $cPValue == 0.090 ] ; then
	pvalueColumn=3
    elif [ $cPValue == 0.080 ] ; then
	pvalueColumn=4
    elif [ $cPValue == 0.070 ] ; then
	pvalueColumn=5
    elif [ $cPValue == 0.060 ] ; then
	pvalueColumn=6
    elif [ $cPValue == 0.050 ] ; then
	pvalueColumn=7
    elif [ $cPValue == 0.040 ] ; then
	pvalueColumn=8
    elif [ $cPValue == 0.030 ] ; then
	pvalueColumn=9
    elif [ $cPValue == 0.020 ] ; then
	pvalueColumn=10
    elif [ $cPValue == 0.010 ] ; then
	pvalueColumn=11
    fi

    nVoxels=$( cat $clustsimPrefix.NN$nn.1D | sed '/^#/d' | grep "^ $vPValue" | awk "{print \$$pvalueColumn}" )

    echo $nVoxels
}

## nearest neighbour 1=touching at faces, 2=faces and edges 3=faces,
## edges and corners, just like in the afni clusterize window
NN=1
## the corrected p-value for the clusters. NB: the trailing zero at the
## end of this p-value is really important, do not omit it
##cPvalue=0.050
case $NN in
    1)
	rmm=1.01
	;;
    2)
	rmm=1.44
	;;
    3)
	rmm=1.75
	;;
esac


#export OMP_NUM_THREADS=8
usedFwhm=4.2

###################if [ 1 == 0 ] ; then

csvFile=parameters.%auc.csv
echo "contrast,fwhm,fLabel,fValueBrikId,fThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,bucketFile" > $csvFile

for contrast in $contrasts ; do

    bucketFilename="$GROUP_DATA/$task.bucket.$groups.$contrast.%auc.masked+tlrc.HEAD"
    ##                         ("pine.bucket", groups, contrast, "%auc.masked+tlrc.HEAD"
    if [ ! -f $bucketFilename ] ; then
	pwd
	echo "*** ERROR cannot find $bucketFilename. Exiting"
	exit
    fi
    
	##fwhm=( $( 3dTstat -mean -prefix stdout: riskygains.$groups.$phase.allSubjects.fwhm.ext.1D\' 2> /dev/null ) )
	##usedFwhm=${fwhm[3]
    
    if [ ! -f CStemp.fwhm${usedFwhm}.NN1.1D ] ; then
	echo "*** Running 3dClustSim"
	#export OMP_NUM_THREADS=10
	if [ -f $GROUP_RESULTS/mask+tlrc.HEAD ] ; then 
	    3dClustSim -mask mask+tlrc.HEAD -fwhm ${usedFwhm} -niml -prefix CStemp.fwhm${usedFwhm}
	else 
	    3dClustSim -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz  -fwhm ${usedFwhm} -niml -prefix CStemp.fwhm${usedFwhm}
	fi
    fi
    
    latestLmeBucketFile=$( pickLatestBucketFile $contrast ) 
    echo "** Most recent bucket file is $latestLmeBucketFile"

    3drefit \
	-atrstring AFNI_CLUSTSIM_NN1 file:CStemp.fwhm${usedFwhm}.NN1.niml \
	-atrstring AFNI_CLUSTSIM_MASK file:CStemp.fwhm${usedFwhm}.mask \
	-atrstring AFNI_CLUSTSIM_NN2 file:CStemp.fwhm${usedFwhm}.NN2.niml \
	-atrstring AFNI_CLUSTSIM_NN3 file:CStemp.fwhm${usedFwhm}.NN3.niml \
	$latestLmeBucketFile
	    
    fValueBrikLabels=$( 3dinfo -label $latestLmeBucketFile | tr "|" "\n" | grep "F-value" 2> /dev/null )
    for fLabel in $fValueBrikLabels ; do
	fValueBrikId=$( 3dinfo -label2index $fLabel $latestLmeBucketFile 2> /dev/null ) 
	echo "### Parameters are:"
	echo "### $latestLmeBucketFile: $fLabel $fValueBrikId"
	## pick the p-values here 
	pValue=0.05
	cPvalue=0.050
	
	nVoxels=$( extractNVoxels $NN $cPvalue $pValue CStemp.fwhm${usedFwhm} ) 
	df=$( extractFStatpars $latestLmeBucketFile $fValueBrikId )
		
	fThreshold=$( cdf -p2t fift $pValue $df | sed 's/t = //' )
	echo "### fwhm = ${usedFwhm}" # fwhm[3]}"
	echo "### fLabel = $fLabel"
	echo "### fValueBrikId = $fValueBrikId"
	echo "### fThreshold = $fThreshold"
	echo "### rmm = $rmm"
	echo "### nVoxels = $nVoxels"
	echo "### df = $df"
	echo "### voxelwise pValue = $pValue"
	echo "### corrected  pValue = $cPvalue"
	
	suffix="fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel.%auc"
	3dclust -1Dformat -savemask clorder.$suffix -nosum -1dindex $fValueBrikId -1tindex $fValueBrikId -1noneg -2thresh -$fThreshold $fThreshold \
	    -dxyz=1 $rmm $nVoxels $latestLmeBucketFile  > clust.$suffix.txt
	if [ -f clorder.$suffix+tlrc.HEAD ] ; then 
	    nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )
	    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename > roiStats.$suffix.txt
	    # 3dTstat -mean -mask clorder.$suffix+tlrc.HEAD -prefix clust.mean.$suffix $bucketFilename 
	    # 3drefit \
	    # 	-atrstring AFNI_CLUSTSIM_NN1 file:CStemp.fwhm${usedFwhm}.NN1.niml \
	    # 	-atrstring AFNI_CLUSTSIM_MASK file:CStemp.fwhm${usedFwhm}.mask \
	    # 	-atrstring AFNI_CLUSTSIM_NN2 file:CStemp.fwhm${usedFwhm}.NN2.niml \
	    # 	-atrstring AFNI_CLUSTSIM_NN3 file:CStemp.fwhm${usedFwhm}.NN3.niml \
	    # 	clust.mean.$suffix+tlrc.HEAD
	    3drefit \
		-atrstring AFNI_CLUSTSIM_NN1 file:CStemp.fwhm${usedFwhm}.NN1.niml \
		-atrstring AFNI_CLUSTSIM_MASK file:CStemp.fwhm${usedFwhm}.mask \
		-atrstring AFNI_CLUSTSIM_NN2 file:CStemp.fwhm${usedFwhm}.NN2.niml \
		-atrstring AFNI_CLUSTSIM_NN3 file:CStemp.fwhm${usedFwhm}.NN3.niml \
		clorder.$suffix+tlrc.HEAD
	else
	    nClusters=0
	    echo "*** WARNING No clusters found!"
	fi
	echo "$contrast,${usedFwhm},$fLabel,$fValueBrikId,$fThreshold,$rmm,$nVoxels,$df,$pValue,$cPvalue,$nClusters,$latestLmeBucketFile" >> $csvFile
    done ## end of for fLabel in $fValueBrikLabels ; do
done ## end of for contrast in $contrasts ; do


## contrasts="emotiveVsNeutral allRemVsAllNotrem"

# echo "*** Now extracting the percent change from the ROI surviving the whole-brain analysis"

# ## now that the results of the GLT LME are clustered, we need to
# ## extract the values in the percent change bucket files for graphing
# # contrast=emotiveVsNeutral

# for contrast in $contrasts ; do 
#     suffix=fwhm${usedFwhm}.$task.$groups.$contrast.group.F-value
#     if [ -f clorder.$suffix+tlrc.HEAD ] ; then

#     # "fearfulVsHappy"            = c("fearful", "happy"),
#     # "fearfulVsNeutral"          = c("fearful", "Neutral"),
#     # "fearfulVsSad"              = c("Fearful", "sad"),
#     # "happyVsNeutral"            = c("Happy", "neutral"),
#     # "happyVsSad"                = c("Happy", "sad"),
#     # "neutralVsSad"              = c("neutral", "sad"),
#     # "allEmotiveVsNeutral"       = c("fearful", "happy", "sad", "neutral"),
#     # "happyRemVsHappyNotrem"     = c("happyRem", "happyNotrem"),
#     # "fearfulRemVsFearfulNotrem" = c("fearfulRem", "fearfulNotrem"),
#     # "neutralRemVsNeutralNotrem" = c("neutralRem", "neutralNotrem"),
#     # "sadRemVsSadNotrem"         = c("sadRem", "sadNotRem"),
#     # "allRemVsAllNotrem"         = c("allRem", "allNotrem")

# 	if [ $contrast == "fearfulVsHappy" ] ; then 
# 	    stimuli="fearful happy"
# 	elif [ $contrast == "fearfulVsNeutral" ] ; then 
# 	    stimuli="fearful Neutral"
# 	elif [ $contrast == "fearfulVsSad" ] ; then 
# 	    stimuli="fearful sad"
# 	elif [ $contrast == "happyVsNeutral" ] ; then 
# 	    stimuli="happy neutral"
# 	elif [ $contrast == "happyVsSad" ] ; then 
# 	    stimuli="happy sad"
# 	elif [ $contrast == "neutralVsSad" ] ; then 
# 	    stimuli="neutral sad"
# 	elif [ $contrast == "allEmotiveVsNeutral" ] ; then 
# 	    stimuli="happy fearful sad neutral"
# 	elif [ $contrast == "happyRemVsHappyNotrem" ] ; then 
# 	    stimuli="happyRem happyNotrem"
# 	elif [ $contrast == "fearfulRemVsFearfulNotrem" ] ; then 
# 	    stimuli="fearfulRem fearfulNotrem"
# 	elif [ $contrast == "neutralRemVsNeutralNotrem" ] ; then 
# 	    stimuli="neutralRem neutralNotrem"
# 	elif [ $contrast == "sadRemVsSadNotrem" ] ; then 
# 	    stimuli="sadRem sadNotRem"
# 	elif [ $contrast == "allRemVsAllNotrem" ] ; then 
# 	    stimuli="happyRem happyNotrem fearfulRem fearfulNotrem sadRem sadNotrem neutralRem neutralNotrem"
# 	else
# 	    exit
# 	fi

# 	for stimulus in $stimuli ; do
# 	    echo "*** Found stimulus $stimulus in the contrast $contrast. Extracting %CS from ROIs"
# 	    3dROIstats  -mask clorder.$suffix+tlrc.HEAD $GROUP_DATA/${task}.bucket.$groups.$stimulus.%cs.REML.masked+tlrc > \
# 		roiStats.pctChange.${stimulus}.stimulus.${suffix}.contrast.$groups.txt
# 	done
#     fi
# done
 
cd $scriptsDir
echo "*** Making cluster location tables"
./cluster2Table.pl --space=mni --force $GROUP_RESULTS

