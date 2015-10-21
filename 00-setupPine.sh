#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM

#set -x 

taskDir=cPine
task=Pine

dataRoot="/Volumes/PROMISEPEGASUS/yangdata/"
scriptsDir="/Volumes/PROMISEPEGASUS/yangdata/"

if [ $# -gt 0 ] ; then
    subjects="$*"
else
    subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.subjectList.txt )"
fi

function reconstructAnatomy {
    subject=$1
    dicomTask=FSPGR_SAG_TI550

    subjectDicomContainerDir=$dataRoot/$subject
    if [ -d $subjectDicomContainerDir ] ; then 
	if [ -f $subjectDicomContainerDir/sinfo.txt ] ; then 
	    sdir=$( cat $subjectDicomContainerDir/sinfo.txt | grep $dicomTask |awk '{print $1}' )
	else 
	    (cd $subjectDicomContainerDir; ../scheck )
	    #echo "*** Could not find $subjectDicomContainerDir/sinfo.txt from which to determine the correct s-directory containing the DICOMS for the $task task"
	    #echo "*** Please provide the s-directory name (enter Ctrl-C to quit at this point if you do not know the correct s-directory)"
	    #read sdir
	    ##echo "** There are ${#sdir} characters in the entered text"
	    ##exit
	fi
	if [ ${#sdir} -gt 0 ] && [ -d  $subjectDicomContainerDir/$sdir ] ; then 
	    echo "*** Resequencing the i files in $subjectDicomContainerDir/$sdir"
	    ( cd $subjectDicomContainerDir/; /Volumes/opt/local/bin/perl $scriptsDir/imseq $sdir )
	    
	    session=$dataRoot/${subject}BRIKS
	    if [ ! -d $session ] ; then 
		echo "*** $session does not exist. Creating it now"
		mkdir -p $session
	    fi

	    prefix="$subject"
	    if [ ! -f $session/$prefix+orig.HEAD ] ; then 
		echo "*** Now creating AFNI HEAD/BRIK for $subject"
		echo "*** Prefix will be $prefix"
		( cd $subjectDicomContainerDir/$sdir ;  to3d -anat -prefix $prefix -session $session "i*" )
	    else 
		echo "*** $session/$prefix+orig.HEAD already exists. Skipping creation of anatomy"
	    fi

	else
	    echo "*** The s-directory $sdir does not exist or you entered a empty string for it when prompted. Cannot reconstruct $task task data for subject ${subject}. Skipping."
	fi
    else 
	echo "*** Cannot find $subjectDicomContainerDir"
	echo "*** Skipping"
    fi

}


function prepareForFunctionalReconstruction {
    subject=$1
    subjectDicomContainerDir=$dataRoot/$subject
    if [ -d $subjectDicomContainerDir ] ; then 

	( cd $subjectDicomContainerDir; mv e*/* ./ )
	( cd $subjectDicomContainerDir; $scriptsDir/scheck )
	
    else 
	echo "*** Cannot find $subjectDicomContainerDir"
	echo "*** Skipping"
    fi
}

function reconstructFunctional {
    subject=$1

    dicomTask="PineNew"
    
    slices="40"
    volumes="430"
    tr="2000"
    
    subjectDicomContainerDir=$dataRoot/$subject
    if [ -d $subjectDicomContainerDir ] ; then 
	if [ -f $subjectDicomContainerDir/sinfo.txt ] ; then 

	    ## if there are multiple matches tha tail -1 command will
	    ## eliminate all but the last one
	    sdir=$( cat $subjectDicomContainerDir/sinfo.txt | grep $dicomTask | tail -1 | awk '{print $1}' )
	else 
	    (cd $subjectDicomContainerDir; ../scheck )
	    #echo "*** Could not find $subjectDicomContainerDir/sinfo.txt from which to determine the correct s-directory containing the DICOMS for the $task task"
	    #echo "*** Please provide the s-directory name (enter Ctrl-C to quit at this point if you do not know the correct s-directory)"
	    #read sdir
	    ##echo "** There are ${#sdir} characters in the entered text"
	    ##exit
	fi
	if [ ${#sdir} -gt 0 ] && [ -d  $subjectDicomContainerDir/$sdir ] ; then 
	    echo "*** Resequencing the i files in $subjectDicomContainerDir/$sdir"
	    ( cd $subjectDicomContainerDir/; /Volumes/opt/local/bin/perl $scriptsDir/imseq $sdir )
	    
	    session=$dataRoot/${subject}BRIKS
	    if [ ! -d $session ] ; then 
		echo "*** $session does not exist. Creating it now"
		mkdir -p $session
	    fi

	    ##prefix="$subject.$task"
	    prefix="$subject$task"
	    echo "*** Now creating AFNI HEAD/BRIK for $subject"
	    echo "*** Prefix will be $prefix"
	    ( cd $subjectDicomContainerDir/$sdir ;  to3d -epan -prefix $prefix$nifti -session $session -time:zt $slices $volumes $tr alt+z "i*" )

	else
	    echo "*** The s-directory $sdir does not exist or you entered a empty string for it when prompted. Cannot reconstruct $task task data for subject ${subject}. Skipping."
	fi
    else 
	echo "*** Cannot find $subjectDicomContainerDir"
	echo "*** Skipping"
    fi


}

function linkAnatomy {
    subject=$1
    echo "*** Linking anatomy for $subject"

    if [ -f $dataRoot/${subject}BRIKS/${subject}+orig.HEAD ] ; then 
	if [ ! -f $dataRoot/${taskDir}/data/${subject}/${subject}+orig.HEAD ] ; then 
	    ln -sf $dataRoot/${subject}BRIKS/${subject}+orig.HEAD $dataRoot/${taskDir}/data/${subject}/
	fi
    fi

    if [ -f $dataRoot/${subject}BRIKS/${subject}+orig.BRIK ] ; then 
	if [ ! -f $dataRoot/${taskDir}/data/${subject}/${subject}+orig.BRIK ] ; then 
	    ln -sf $dataRoot/${subject}BRIKS/${subject}+orig.BRIK $dataRoot/${taskDir}/daeta/${subject}/
	fi
    else 
	if [ -f $dataRoot/${subject}BRIKS/${subject}+orig.BRIK.gz ] ; then 
		## assume that the gzipped version exists
	    if [ ! -f $dataRoot/${taskDir}/data/${subject}/${subject}+orig.BRIK.gz ] ; then 
		ln -sf $dataRoot/${subject}BRIKS/${subject}+orig.BRIK.gz $dataRoot/${taskDir}/data/${subject}/
	    fi
	fi
    fi
}

function linkFunctional {
    subject=$1
    echo "*** Linking functional for $subject"

    if [ -f  $dataRoot/${subject}BRIKS/${subject}${task}+orig.HEAD ] ; then 
	if [ ! -f $dataRoot/${taskDir}/data/${subject}/${subject}${task}+orig.HEAD ] ; then 
	    ln -sf $dataRoot/${subject}BRIKS/${subject}${task}+orig.HEAD $dataRoot/${taskDir}/data/${subject}/
	fi
    fi

    if [ -f $dataRoot/${subject}BRIKS/${subject}${task}+orig.BRIK ] ; then 
	if [ ! -f $dataRoot/${taskDir}/data/${subject}/${subject}${task}+orig.BRIK ] ; then 
	    ln -sf $dataRoot/${subject}BRIKS/${subject}${task}+orig.BRIK $dataRoot/${taskDir}/data/${subject}/
	fi
    else 
	if [ -f $dataRoot/${subject}BRIKS/${subject}${task}+orig.BRIK.gz ] ; then 
		## assume that the gzipped version exists
	    if [ ! -f $dataRoot/${taskDir}/data/${subject}/ ] ; then 
		ln -sf $dataRoot/${subject}BRIKS/${subject}${task}+orig.BRIK.gz $dataRoot/${taskDir}/data/${subject}/
	    fi
	fi
    fi
}


noDataDir=""
existingDataDirs=""

for subjectNumber in $subjects ; do
    for timepoint in A ; do
	subjectNumber="${subjectNumber%_*}_${timepoint}"

	echo "####################################################################################################"
	echo "### Timepoint $timepoint: $subjectNumber"
	
	if [ -d $dataRoot/$subjectNumber ] ; then 
	    existingDataDirs="$existingDataDirs $subjectNumber"
	    if [ ! -d $dataRoot/${taskDir}/data/$subjectNumber ] ; then 
		mkdir $dataRoot/${taskDir}/data/$subjectNumber
	    fi
	    
	    if [ ! -f $dataRoot/$subjectNumber/sinfo.txt ] ; then 
		echo "*** Preparing for reconstruction of functional for $subjectNumber"
		prepareForFunctionalReconstruction  $subjectNumber
	    fi
	    
	    if [ -f $dataRoot/$subjectNumber+orig.HEAD ] ; then 
		linkAnatomy $subjectNumber
	    else
		echo "*** Trying to reconstruct anatomy for $subjectNumber"
		reconstructAnatomy $subjectNumber
		linkAnatomy $subjectNumber
	    fi
	else 
	    noDataDir="$noDataDir $subjectNumber"
	    echo "*** $dataRoot/$subjectNumber does not exist. Skipping"
	fi
	
	### Now deal with the functional 

	#rm -f $dataRoot/${subjectNumber}BRIKS/${subjectNumber}${task}+orig.*

     	if [ -f $dataRoot/${subjectNumber}BRIKS/${subjectNumber}${task}+orig.HEAD ] ; then 
	    linkFunctional $subjectNumber
	else
	    ## no functional? try to do reconstruction
	    
	    echo "*** Trying to reconstruct functional for $subjectNumber"
	    reconstructFunctional $subjectNumber
	    linkFunctional $subjectNumber
	fi
	
    done ## end of for timepoint in C D ; do
done ## end of for subjectNumber in $subjects ; do

echo "####################################################################################################"
echo "Found data directories for:"
echo $existingDataDirs
echo "####################################################################################################"
echo "No data directories for:"
echo $noDataDir

