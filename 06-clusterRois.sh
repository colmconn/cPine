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

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad" # allEmotiveVsNeutral happyRemVsHappyNotrem fearfulRemVsFearfulNotrem neutralRemVsNeutralNotrem sadRemVsSadNotrem allRemVsAllNotrem"

cd $GROUP_RESULTS
groups="mddAndCtrl"

function pickLatestBucketFile {
    local contrast=$1
    ##restingstate.$groups.$contrast.REML.lme.bucket.*+tlrc.HEAD

    latest=$( ls -1t $task.lme.bucket.$groups.$contrast.REML.*+tlrc.HEAD | grep -v masked | head -1 ) 
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

GETOPT_OPTIONS=$( $GETOPT  -o "is:n:" --longoptions "useInherentSmoothness,svc:,nn:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

doSvc=0
useInherentSmoothness=0
# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--svc)
	    doSvc=1; 
	    ## masks="HarvardOxford-sub-maxprob-thr25-3mm-bilateral-amygdala+tlrc.HEAD sgacc.bilateral.3mm+tlrc.HEAD"
	    ## masks="bilateralAmygdalaeAndSgaccRois+tlrc.HEAD"
	    masks="$2"
	    echo "*** Doing small volume correction for the following masks: $masks"
	    shift 2 ;;
	-n|--nn)
	    NN=$2; 
	    shift 2 ;;
	-i|--useInherentSmoothness ) 
	    useInherentSmoothness=1; 
	    shift ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ "x$NN" == "x" ] ; then 
    ## nearest neighbour 1=touching at faces, 2=faces and edges 3=faces,
    ## edges and corners, just like in the afni clusterize window
    NN=1
fi
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

    *) 
	echo "Unknown value ($NN) for NN. Exiting."
	exit 2 ;;
esac


#export OMP_NUM_THREADS=8

if [ $useInherentSmoothness -eq 1 ] ; then 
    if [ -f $GROUP_DATA/pine.ctrlOnly.REML.fwhmEstimates.tab ] ; then 
	usedFwhm=$( cat $GROUP_DATA/pine.ctrlOnly.REML.fwhmEstimates.tab | gawk '{print $4}' | gawk '{s+=$1}END{print s/NR}' )
    else 
	usedFwhm=4.2
    fi
else 
    usedFwhm=4.2
fi

echo "*** Correcting for $usedFwhm mm smoothness"

###################if [ 1 == 0 ] ; then

if  [ $doSvc -eq 1 ] ; then 
    csvFile=parameters.svc.csv
    echo "contrast,fwhm,fLabel,fValueBrikId,fThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,mask,bucketFile" > $csvFile
else 
    csvFile=parameters.csv
    echo "contrast,fwhm,fLabel,fValueBrikId,fThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,mask,bucketFile" > $csvFile
fi

for contrast in $contrasts ; do

    echo "####################################################################################################"
    echo "### Contrast is: $contrast"

    bucketFilename="$GROUP_DATA/$task.bucket.$groups.$contrast.REML.masked+tlrc.HEAD"
    ctrlOnlyBucketFilename="$GROUP_DATA/$task.bucket.ctrlOnly.$contrast.REML.masked+tlrc.HEAD"
    mddOnlyBucketFilename="$GROUP_DATA/$task.bucket.mddOnly.$contrast.REML.masked+tlrc.HEAD"

    if [ ! -f $bucketFilename ] ; then
	pwd
	echo "*** ERROR cannot find $bucketFilename. Exiting"
	exit
    fi
    
    ## fwhm=( $( 3dTstat -mean -prefix stdout: riskygains.$groups.$phase.allSubjects.fwhm.ext.1D\' 2> /dev/null ) )
    ## usedFwhm=${fwhm[3]
  
    ## delete all the masked versions of the latest LME bucket file to
    ## prevent the pickLatestBucketFile function from finding them instead
    ## of the original non-masked version
    ## rm -f $task.lme.bucket.$groups.$contrast.REML.*+tlrc.HEAD

    latestLmeBucketFile=$( pickLatestBucketFile $contrast ) 
    echo "### Most recent bucket file is $latestLmeBucketFile"
    fValueBrikLabels=$( 3dinfo -label $latestLmeBucketFile | tr "|" "\n" | grep "F-value" | grep -i "group" 2> /dev/null )

    if  [ $doSvc -eq 1 ] ; then 

	## ####################################################################################################
	## A-priori small volume masks
	## ####################################################################################################

	#csvFile=parameters.svc.csv

	for mask in $masks ; do
	    echo "*** Mask is: $mask"
	    maskName=${mask%%+*}

	    cstempPrefix=CStemp.fwhm${usedFwhm}.${maskName}
	    if [ ! -f ${cstempPrefix}.NN1.1D ] ; then
		echo "*** Running 3dClustSim"
		#export OMP_NUM_THREADS=10
		3dClustSim -mask $ROIS_DIR/$mask -fwhm ${usedFwhm} -niml -prefix ${cstempPrefix}
	    else 
		echo "Found $cstempPrefix.NN$NN.1D. Not running 3dClustSim"
	    fi

	    for fLabel in $fValueBrikLabels ; do
		fValueBrikId=$( 3dinfo -label2index $fLabel $latestLmeBucketFile 2> /dev/null ) 
		echo "### Parameters are:"
		echo "### $latestLmeBucketFile: $fLabel $fValueBrikId"
		## pick the p-values here 
		pValue=0.05
		cPvalue=0.050
	
		nVoxels=$( extractNVoxels $NN $cPvalue $pValue ${cstempPrefix} ) 
		df=$( extractFStatpars $latestLmeBucketFile $fValueBrikId )
		
		fThreshold=$( cdf -p2t fift $pValue $df | sed 's/t = //' )
		##echo "### fwhm = ${fwhm[3]}"
		echo "### fwhm = ${fwhm}"
		echo "### fLabel = $fLabel"
		echo "### fValueBrikId = $fValueBrikId"
		echo "### fThreshold = $fThreshold"
		echo "### rmm = $rmm"
		echo "### nVoxels = $nVoxels"
		echo "### df = $df"
		echo "### voxelwise pValue = $pValue"
		echo "### corrected  pValue = $cPvalue"
	
		suffix=fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel.${maskName}
		mddOnlySuffix=fwhm${usedFwhm}.$task.mddOnly.$contrast.$fLabel.${maskName}
		ctrlOnlySuffix=fwhm${usedFwhm}.$task.ctrlOnly.$contrast.$fLabel.${maskName}
		
		## mask the bucket file so we don't get clusters outside the regions of interest
		
		# echo
		# echo
		# echo
		# echo "mask: ${mask}"
		# echo "Latest lme bucket file: ${latestLmeBucketFile}"
		# echo "Latest LME bucket file: ${latestLmeBucketFile%%+*}"
		maskedLatestLmeBucketFilePrefix=${latestLmeBucketFile%%+*}.${mask%%+*}.masked
		# echo "Latest masked LME bucket file: ${maskedLatestLmeBucketFilePrefix}"
		
		3dcalc -datum float -a $latestLmeBucketFile -b $ROIS_DIR/$mask -expr "a*b" -prefix $maskedLatestLmeBucketFilePrefix
		
		3dclust -1Dformat -savemask clorder.$suffix -nosum -1dindex $fValueBrikId -1tindex $fValueBrikId -1noneg -2thresh -$fThreshold $fThreshold \
		    -dxyz=1 $rmm $nVoxels ${maskedLatestLmeBucketFilePrefix}+tlrc.  > clust.$suffix.txt
		##rm -f ${maskedLatestLmeBucketFilePrefix}*
		if [ -f clorder.$suffix+tlrc.HEAD ] ; then 
		    nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )
		    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename         > roiStats.$suffix.txt
		    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $ctrlOnlyBucketFilename > roiStats.$ctrlOnlySuffix.txt
		    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $mddOnlyBucketFilename  > roiStats.$mddOnlySuffix.txt
		    
		    3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${maskedLatestLmeBucketFilePrefix}+tlrc.HEAD\[$fValueBrikId\]         > roiStats.$suffix.meanFValue.txt

		    # 3dTstat -mean -mask clorder.$suffix+tlrc.HEAD -prefix clust.mean.$suffix $bucketFilename 
		    # 3drefit \
		    # 	-atrstring AFNI_CLUSTSIM_NN1 file:CStemp.fwhm${usedFwhm}.NN1.niml \
		    # 	-atrstring AFNI_CLUSTSIM_MASK file:CStemp.fwhm${usedFwhm}.mask \
		    # 	-atrstring AFNI_CLUSTSIM_NN2 file:CStemp.fwhm${usedFwhm}.NN2.niml \
		    # 	-atrstring AFNI_CLUSTSIM_NN3 file:CStemp.fwhm${usedFwhm}.NN3.niml \
		    # 	clust.mean.$suffix+tlrc.HEAD

		    3drefit \
			-atrstring AFNI_CLUSTSIM_NN1  file:${cstempPrefix}.NN1.niml \
			-atrstring AFNI_CLUSTSIM_MASK file:${cstempPrefix}.mask \
			-atrstring AFNI_CLUSTSIM_NN2  file:${cstempPrefix}.NN2.niml \
			-atrstring AFNI_CLUSTSIM_NN3  file:${cstempPrefix}.NN3.niml \
			clorder.$suffix+tlrc.HEAD
		    3drefit -cmap INT_CMAP clorder.$suffix+tlrc.HEAD
		else
		    nClusters=0
		    echo "*** WARNING No clusters found!"
		fi
		echo "$contrast,${usedFwhm},$fLabel,$fValueBrikId,$fThreshold,$rmm,$nVoxels,$df,$pValue,$cPvalue,$nClusters,$mask,$latestLmeBucketFile" >> $csvFile
	    done ## end of for fLabel in $fValueBrikLabels ; do
	    
	done ## end of for mask in ...
    else 
	## ####################################################################################################
	## Whole brain mask
	## ####################################################################################################

	#csvFile=parameters.csv
	#echo "contrast,fwhm,fLabel,fValueBrikId,fThreshold,rmm,nVoxels,df,pValue,cPvalue,nClusters,mask,bucketFile" > $csvFile

	cstempPrefix=CStemp.fwhm${usedFwhm}
	if [ ! -f ${cstempPrefix}.NN1.1D ] ; then
	    echo "*** Running 3dClustSim"
	    #export OMP_NUM_THREADS=10
	    if [ -f $GROUP_RESULTS/mask+tlrc.HEAD ] ; then 
		3dClustSim -mask mask+tlrc.HEAD -fwhm ${usedFwhm} -niml -prefix ${cstempPrefix}
	    else 
		3dClustSim -mask $MDD_STANDARD/MNI152_T1_3mm_brain_mask.nii.gz  -fwhm ${usedFwhm} -niml -prefix ${cstempPrefix}
	    fi
	fi
	3drefit \
	    -atrstring AFNI_CLUSTSIM_NN1  file:${cstempPrefix}.NN1.niml \
	    -atrstring AFNI_CLUSTSIM_MASK file:${cstempPrefix}.mask \
	    -atrstring AFNI_CLUSTSIM_NN2  file:${cstempPrefix}.NN2.niml \
	    -atrstring AFNI_CLUSTSIM_NN3  file:${cstempPrefix}.NN3.niml \
	    $latestLmeBucketFile
    
	fValueBrikLabels=$( 3dinfo -label $latestLmeBucketFile | tr "|" "\n" | grep "F-value" | grep -i "group" 2> /dev/null )
	for fLabel in $fValueBrikLabels ; do
	    fValueBrikId=$( 3dinfo -label2index $fLabel $latestLmeBucketFile 2> /dev/null ) 
	    echo "### Parameters are:"
	    echo "### $latestLmeBucketFile: $fLabel $fValueBrikId"
	    ## pick the p-values here 
	    pValue=0.05
	    cPvalue=0.050
	
	    nVoxels=$( extractNVoxels $NN $cPvalue $pValue ${cstempPrefix} ) 
	    df=$( extractFStatpars $latestLmeBucketFile $fValueBrikId )
		
	    fThreshold=$( cdf -p2t fift $pValue $df | sed 's/t = //' )
	    echo "### fwhm = ${fwhm[3]}"
	    echo "### fLabel = $fLabel"
	    echo "### fValueBrikId = $fValueBrikId"
	    echo "### fThreshold = $fThreshold"
	    echo "### rmm = $rmm"
	    echo "### nVoxels = $nVoxels"
	    echo "### df = $df"
	    echo "### voxelwise pValue = $pValue"
	    echo "### corrected  pValue = $cPvalue"
	
	    suffix=fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel
	    mddOnlySuffix=fwhm${usedFwhm}.$task.mddOnly.$contrast.$fLabel
	    ctrlOnlySuffix=fwhm${usedFwhm}.$task.ctrlOnly.$contrast.$fLabel
	    
	    3dclust -1Dformat -savemask clorder.$suffix -nosum -1dindex $fValueBrikId -1tindex $fValueBrikId -1noneg -2thresh -$fThreshold $fThreshold \
		-dxyz=1 $rmm $nVoxels $latestLmeBucketFile  > clust.$suffix.txt
	    if [ -f clorder.$suffix+tlrc.HEAD ] ; then 
		nClusters=$( 3dBrickStat -max clorder.$suffix+tlrc.HEAD 2> /dev/null | tr -d ' ' )
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $bucketFilename         > roiStats.$suffix.txt
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $ctrlOnlyBucketFilename > roiStats.$ctrlOnlySuffix.txt
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD $mddOnlyBucketFilename  > roiStats.$mddOnlySuffix.txt
	    
		3dROIstats -nobriklab -mask clorder.$suffix+tlrc.HEAD ${latestLmeBucketFile}\[$fValueBrikId\]         > roiStats.$suffix.averageFvalue.txt

		# 3dTstat -mean -mask clorder.$suffix+tlrc.HEAD -prefix clust.mean.$suffix $bucketFilename 
		# 3drefit \
		# 	-atrstring AFNI_CLUSTSIM_NN1 file:CStemp.fwhm${usedFwhm}.NN1.niml \
		# 	-atrstring AFNI_CLUSTSIM_MASK file:CStemp.fwhm${usedFwhm}.mask \
		# 	-atrstring AFNI_CLUSTSIM_NN2 file:CStemp.fwhm${usedFwhm}.NN2.niml \
		# 	-atrstring AFNI_CLUSTSIM_NN3 file:CStemp.fwhm${usedFwhm}.NN3.niml \
		# 	clust.mean.$suffix+tlrc.HEAD
		3drefit \
		    -atrstring AFNI_CLUSTSIM_NN1  file:${cstempPrefix}.NN1.niml \
		    -atrstring AFNI_CLUSTSIM_MASK file:${cstempPrefix}.mask \
		    -atrstring AFNI_CLUSTSIM_NN2  file:${cstempPrefix}.NN2.niml \
		    -atrstring AFNI_CLUSTSIM_NN3  file:${cstempPrefix}.NN3.niml \
		    clorder.$suffix+tlrc.HEAD
		3drefit -cmap INT_CMAP clorder.$suffix+tlrc.HEAD
	    else
		nClusters=0
		echo "*** WARNING No clusters found!"
	    fi
	    echo "$contrast,${usedFwhm},$fLabel,$fValueBrikId,$fThreshold,$rmm,$nVoxels,$df,$pValue,$cPvalue,$nClusters,$mask,$latestLmeBucketFile" >> $csvFile
	done ## end of for fLabel in $fValueBrikLabels ; do
    fi ## end of if [ $doSvc ] ; then
done ## end of for contrast in $contrasts ; do


## contrasts="emotiveVsNeutral allRemVsAllNotrem"

echo "*** Now extracting the percent change from the ROI surviving the whole-brain analysis"

## now that the results of the GLT LME are clustered, we need to
## extract the values in the percent change bucket files for graphing
# contrast=emotiveVsNeutral

for contrast in $contrasts ; do 
    if [ $contrast == "fearfulVsHappy" ] ; then 
	stimuli="fearful happy"
    elif [ $contrast == "fearfulVsNeutral" ] ; then 
	stimuli="fearful neutral"
    elif [ $contrast == "fearfulVsSad" ] ; then 
	stimuli="fearful sad"
    elif [ $contrast == "happyVsNeutral" ] ; then 
	stimuli="happy neutral"
    elif [ $contrast == "happyVsSad" ] ; then 
	stimuli="happy sad"
    elif [ $contrast == "neutralVsSad" ] ; then 
	stimuli="neutral sad"
    elif [ $contrast == "allEmotiveVsNeutral" ] ; then 
	stimuli="happy fearful sad neutral"
    elif [ $contrast == "happyRemVsHappyNotrem" ] ; then 
	stimuli="happyRem happyNotrem"
    elif [ $contrast == "fearfulRemVsFearfulNotrem" ] ; then 
	stimuli="fearfulRem fearfulNotrem"
    elif [ $contrast == "neutralRemVsNeutralNotrem" ] ; then 
	stimuli="neutralRem neutralNotrem"
    elif [ $contrast == "sadRemVsSadNotrem" ] ; then 
	stimuli="sadRem sadNotRem"
    elif [ $contrast == "allRemVsAllNotrem" ] ; then 
	stimuli="happyRem happyNotrem fearfulRem fearfulNotrem sadRem sadNotrem neutralRem neutralNotrem"
    else
	exit
    fi

    if [ $doSvc  -eq 1 ] ; then
	for mask in $masks ; do
	    echo "### Mask is: $mask"
	    maskName=${mask%%+*}
	    
	    suffix=fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel.${maskName}

	    if [ -f clorder.$suffix+tlrc.HEAD ] ; then

		for stimulus in $stimuli ; do
		    echo "*** Found stimulus $stimulus in the contrast $contrast. Extracting %CS from ROIs"
		    3dROIstats  -mask clorder.$suffix+tlrc.HEAD $GROUP_DATA/${task}.bucket.$groups.$stimulus.%cs.REML.masked+tlrc > \
			roiStats.pctChange.${stimulus}.stimulus.${suffix}.contrast.$groups.txt
		done
	    fi
	done
    else 
	suffix=fwhm${usedFwhm}.$task.$groups.$contrast.group.F-value	
	
	if [ -f clorder.$suffix+tlrc.HEAD ] ; then

	    for stimulus in $stimuli ; do
		echo "*** Found stimulus $stimulus in the contrast $contrast. Extracting %CS from ROIs"
		3dROIstats  -mask clorder.$suffix+tlrc.HEAD $GROUP_DATA/${task}.bucket.$groups.$stimulus.%cs.REML.masked+tlrc > \
		    roiStats.pctChange.${stimulus}.stimulus.${suffix}.contrast.$groups.txt
	    done
	fi
    fi
done

cd $scriptsDir
echo "*** Making cluster location tables"
./cluster2Table.pl --space=mni --force $GROUP_RESULTS

