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

CONFIG_DATA=$DATA/config

scriptsDir=${ROOT}/scripts
logDir=${DATA}/log

task=pine

contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad allEmotiveVsNeutral happyRemVsHappyNotrem fearfulRemVsFearfulNotrem neutralRemVsNeutralNotrem sadRemVsSadNotrem allRemVsAllNotrem"

usedFwhm=4.2
fLabel="group.F-value"

cd $GROUP_RESULTS
groups="mddAndCtrl"


GETOPT_OPTIONS=$( $GETOPT  -o "s:" --longoptions "svc:" -n ${programName} -- "$@" )
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

	--) 
	    shift ; break ;;
	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if  [ $doSvc -eq 1 ] ; then 
    for mask in $masks ; do

	echo "*** Mask is: $mask"
	maskName=${mask%%+*}
	    
	for contrast in $contrasts ; do
	    ## clorder.fwhm4.2.pine.mddAndCtrl.fearfulVsNeutral.group.F-value.apriori.rois.noVmpfc.3mm+tlrc.HEAD 

	    suffix=fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel.${maskName}
	    clorderFile=clorder.$suffix+tlrc.HEAD

	    seedConfigFile=$CONFIG_DATA/ppi_seeds/${contrast}.${maskName}.seeds.txt
	    [[ -f $seedConfigFile ]] && rm -f $seedConfigFile

	    if [ -f $clorderFile ] ; then 
		minMax=( $( 3dBrickStat -min -max $clorderFile 2> /dev/null ) )
		
		echo "*** $clorderFile has ${minMax[1]} clusters"
		## clusters are numbered starting at 1.  
		(( minMax[0] = ${minMax[0]} + 1 ))

		[[ ! -d $GROUP_RESULTS/ppi_seeds/$contrast.${maskName} ]] && mkdir -p $GROUP_RESULTS/ppi_seeds/$contrast.${maskName}
		##clusterPrefix=clorder.$task.$groups.$contrast.$fLabel.roi
		clusterPrefix=roi
		for (( ii=${minMax[0]}; ii <= ${minMax[1]}; ii=ii+1 )) ; do
		    3dcalc -a $clorderFile -prefix ppi_seeds/$contrast.${maskName}/${clusterPrefix}${ii} -expr "equals(a, $ii)"
		    
		    echo "\$GroupResultsDir/ppi_seeds/$contrast.${maskName}/${clusterPrefix}${ii}+tlrc.HEAD" >> $seedConfigFile
		done


	    else
		echo "*** No such file $clorderFile. Skipping."
	    fi
	done

    done
else 

    for contrast in $contrasts ; do

	suffix=fwhm${usedFwhm}.$task.$groups.$contrast.$fLabel
	clorderFile=clorder.$suffix+tlrc.HEAD

	seedConfigFile=$CONFIG_DATA/ppi_seeds/${contrast}.seeds.txt
	[[ -f $seedConfigFile ]] && rm -f $seedConfigFile

	if [ -f $clorderFile ] ; then 
	    minMax=( $( 3dBrickStat -min -max $clorderFile 2> /dev/null ) )
	    
	    echo "*** $clorderFile has ${minMax[1]} clusters"
	    ## clusters are numbered starting at 1.  
	    (( minMax[0] = ${minMax[0]} + 1 ))

	    [[ ! -d $GROUP_RESULTS/ppi_seeds/$contrast ]] && mkdir -p $GROUP_RESULTS/ppi_seeds/$contrast
	    ##clusterPrefix=clorder.$task.$groups.$contrast.$fLabel.roi
	    clusterPrefix=roi
	    for (( ii=${minMax[0]}; ii <= ${minMax[1]}; ii=ii+1 )) ; do
		3dcalc -a $clorderFile -prefix ppi_seeds/$contrast/${clusterPrefix}${ii} -expr "equals(a, $ii)"
		
		echo "\$GroupResultsDir/ppi_seeds/$contrast/${clusterPrefix}${ii}+tlrc.HEAD" >> $seedConfigFile
	    done


	else
	    echo "*** No such file $clorderFile. Skipping."
	fi
    done

fi
