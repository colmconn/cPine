#!/bin/bash

## hello world

set -x 

dataRoot="/Volumes/PROMISEPEGASUS/yangdata/"
scriptsDir="/Volumes/PROMISEPEGASUS/yangdata/"
dicomTask="PineNew"
task="Pine"

slices="40"
volumes="430"
tr="2000"

## uncomment the next line to get a NIfTI file instead of HEAD/BRIK"
##nifti=".nii"

## the next variable is a space seperated list of the subjects to be reconstructed
## and example list of subjects woudl look like
##subjectList="332_A 111_A 222_A 333_A"

#subjectlist="332_A" # 340_A 343_A 356_A 358_A 359_A 360_A 361_A 362_A 366_A 339_A 349_A"
##subjectList="340_A 343_A 356_A 358_A 359_A 360_A 361_A 362_A 366_A 339_A 349_A"

subjectList="382_A"

## edit the next line add your subjects
##subjectList="365_A"

for subject in $subjectList ; do

    subjectDicomContainerDir=$dataRoot/$subject
    if [ -d $subjectDicomContainerDir ] ; then 
	if [ -f $subjectDicomContainerDir/sinfo.txt ] ; then 
	    sdir=$( cat $subjectDicomContainerDir/sinfo.txt | grep $dicomTask |awk '{print $1}' )
	else 
	    (cd $subjectDicomContainerDir; ../scheck )
	    sdir=$( cat $subjectDicomContainerDir/sinfo.txt | grep $dicomTask |awk '{print $1}' )
	    #echo "*** Could not find $subjectDicomContainerDir/sinfo.txt from which to determine the correct s-directory containing the DICOMS for the $task task"
	    #echo "*** Please provide the s-directory name (enter Ctrl-C to quit at this point if you do not know the correct s-directory)"
	    #read sdir
	    ##echo "** There are ${#sdir} characters in the entered text"
	    ##exit
	fi
	if [ ${#sdir} -gt 0 ] && [ -d  $subjectDicomContainerDir/$sdir ] ; then 
	    echo "*** Resequencing the i files in $subjectDicomContainerDir/$sdir"
	    ( cd $subjectDicomContainerDir/; $scriptsDir/imseq $sdir )
	    
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

done
