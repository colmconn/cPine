#!/bin/bash

# set -x

ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
scriptsDir=${ROOT}/scripts

#subjects=$( cd /Volumes/PROMISEPEGASUS/yangdata/cPine/data; \ls -1d *_[ACD] )

## problem subjects: 105_A 108_A 114_A 114_C 117_A 117_C 121_A 139_A 144_C 156_A 161_C 301_C 302_A 332_A 336_C 339_C 339_D 349_D 354_A 356_A 357_A 365_C 366_C

## controls only
#ctrl_subjects="105_A 107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 350_A 354_A 357_A 365_A 369_A 374_A 377_A 379_A 380_A"
#ctrl_subjects="107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 354_A 357_A 367_A 369_A 374_A 377_A" 

ctrl_subjects="$( cat ../data/config/control.subjectList.txt )"

## MDD subjects only
#mdd_subjects="111_A 112_A 113_A 114_A 117_A 134_A 144_A 147_A 158_A 160_A 161_A 300_A 301_A 304_A 316_A 324_A 330_A 332_A 336_A 340_A 343_A 356_A 358_A 359_A 360_A 361_A 362_A 366_A 323_A 339_A 349_A"
mdd_subjects="111_A 112_A 114_A 132_A 134_A 144_A 147_A 158_A 160_A 161_A 164_A 300_A 304_A 313_A 316_A 317_A 318_A 319_A 320_A 323_A 324_A 330_A 331_A 332_A 336_A 339_A 343_A 349_A 356_A 358_A 359_A 360_A 361_A 362_A 363_A 366_A 371_A 372_A 373_A 376_A"

mdd_subjects="$( cat ../data/config/mdd.subjectList.txt )"
subjects="$ctrl_subjects $mdd_subjects"
#subjects="$*"

cd /Volumes/PROMISEPEGASUS/yangdata/cPine/data/regressors/

for subject in $subjects ; do
    echo "####################################################################################################"
    echo "### Processing subject : $subject"

    $scriptsDir/parsePineData.pl -s=${subject} --regressor=onsetOnly

    regressorFile=/Volumes/PROMISEPEGASUS/yangdata/cPine/data/regressors/${subject}.acquisition.onsetOnly.regressors.tab
    if [ -f $regressorFile ] ; then
	cat $regressorFile | awk '{print $6}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.fearful.onsetOnly.1D
        cat $regressorFile | awk '{print $7}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.fearful.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $8}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.fearful.not.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $9}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.fearful.omission.onsetOnly.1D

        cat $regressorFile | awk '{print $11}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.happy.onsetOnly.1D
        cat $regressorFile | awk '{print $12}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.happy.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $13}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.happy.not.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $14}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.happy.omission.onsetOnly.1D

        cat $regressorFile | awk '{print $16}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.neutral.onsetOnly.1D
        cat $regressorFile | awk '{print $17}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.neutral.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $18}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.neutral.not.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $19}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.neutral.omission.onsetOnly.1D

        cat $regressorFile | awk '{print $20}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.sad.onsetOnly.1D
        cat $regressorFile | awk '{print $21}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.sad.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $22}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.sad.not.remembered.onsetOnly.1D
        cat $regressorFile | awk '{print $23}' | sed 1d | grep -v "9999" | tr '\n' ' ' > $subject.sad.omission.onsetOnly.1D

	for rFile in \
	    $subject.fearful.onsetOnly.1D \
            $subject.fearful.remembered.onsetOnly.1D \
            $subject.fearful.not.remembered.onsetOnly.1D \
            $subject.fearful.omission.onsetOnly.1D \
	    $subject.happy.onsetOnly.1D \
	    $subject.happy.remembered.onsetOnly.1D \
	    $subject.happy.not.remembered.onsetOnly.1D \
	    $subject.happy.omission.onsetOnly.1D \
	    $subject.neutral.onsetOnly.1D \
	    $subject.neutral.remembered.onsetOnly.1D \
	    $subject.neutral.not.remembered.onsetOnly.1D \
	    $subject.neutral.omission.onsetOnly.1D \
	    $subject.sad.onsetOnly.1D \
	    $subject.sad.remembered.onsetOnly.1D \
	    $subject.sad.not.remembered.onsetOnly.1D \
	    $subject.sad.omission.onsetOnly.1D 
	do 
	    #binaryFile=${rFile#*.}
	    #binaryFile=${subject}.acquisition.${binaryFile%.tab}.binary.1D
	    
	    #echo "Converting $rFile to binary"
	    eval $( stat -s $rFile )
            if [[ $st_size -eq 0 ]] ; then 
		echo "*" > $rFile
	    fi
	    #timing_tool.py -timing $rFile -timing_to_1D $binaryFile -tr 2 -min_frac 0.5 -run_len 860
	done
	

    else
	nonExistantSubjects="$nonExistantSubjects $subject"
	#echo "$subject: DOES NOT EXIST"
    fi
done

echo "Regressors for the the following subjects do not exist:"
echo "${nonExistantSubjects:-None}"

nonExistantSubjects=""

## now check that the regressors exist
for subject in $subjects ; do

    $scriptsDir/parsePineData.pl -s=${subject} --regressor=onsetDuration

    regressorFile=/Volumes/PROMISEPEGASUS/yangdata/cPine/data/regressors/${subject}.acquisition.onsetDuration.regressors.tab
    if [ -f $regressorFile ] ; then

    	## These are the column indices of the regressors 
    	#   0     1       InstrA
        #   1     2       InstrB
        #   2     3       InstrC
        #   3     4       InstrD
        #   4     5       InstrR
        #   5     6       Fearful 
        #   6     7       Fearful Remembered
        #   7     8       Fearful Not Remembered
        #   8     9       Fearful Omission
        #   9     10      Fixation
        #   10    11      Happy
        #   11    12      Happy Remembered
        #   12    13      Happy Not Remembered
        #   13    14      Happy Omission
        #   14    15      Instructions
        #   15    16      Neutral
        #   16    17      Neutral Remembered
        #   17    18      Neutral Not Remembered 
        #   18    19      Neutral Omission
        #   19    20      Sad
        #   20    21      Sad Remembered
        #   21    22      Sad Not Remembered 
        #   22    23      Sad Omission

    	## the sed 1d deletes the first line which contains the header
    	## of the column, waver does not like the header

        cat $regressorFile | awk '{print $6}' | sed 1d > rm_$subject.fearful.onsetDuration.tab
        cat $regressorFile | awk '{print $7}' | sed 1d > rm_$subject.fearful.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $8}' | sed 1d > rm_$subject.fearful.not.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $9}' | sed 1d > rm_$subject.fearful.omission.onsetDuration.tab

        cat $regressorFile | awk '{print $11}' | sed 1d > rm_$subject.happy.onsetDuration.tab
        cat $regressorFile | awk '{print $12}' | sed 1d > rm_$subject.happy.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $13}' | sed 1d > rm_$subject.happy.not.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $14}' | sed 1d > rm_$subject.happy.omission.onsetDuration.tab

        cat $regressorFile | awk '{print $16}' | sed 1d > rm_$subject.neutral.onsetDuration.tab
        cat $regressorFile | awk '{print $17}' | sed 1d > rm_$subject.neutral.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $18}' | sed 1d > rm_$subject.neutral.not.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $19}' | sed 1d > rm_$subject.neutral.omission.onsetDuration.tab

        cat $regressorFile | awk '{print $20}' | sed 1d > rm_$subject.sad.onsetDuration.tab
        cat $regressorFile | awk '{print $21}' | sed 1d > rm_$subject.sad.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $22}' | sed 1d > rm_$subject.sad.not.remembered.onsetDuration.tab
        cat $regressorFile | awk '{print $23}' | sed 1d > rm_$subject.sad.omission.onsetDuration.tab
        
        echo "*** Wavering for $subject"

        for g in \
    	    fearful fearful.remembered fearful.not.remembered fearful.omission \
    	    happy happy.remembered happy.not.remembered happy.omission \
    	    neutral neutral.remembered neutral.not.remembered neutral.omission \
    	    sad sad.remembered sad.not.remembered sad.omission ; do

    	    ## 430 is the number of TRs in the scan, so match the length of the regressor to that
            waver -dt 2 -GAM -tstim $( cat rm_$subject.${g}.onsetDuration.tab ) -numout 430 > rm_$subject.${g}.wavered.1D

        done
        
        ## response needs to be conormalised
        
        echo "*** Normalising for $subject"
        # co normalise the response as it can vary in length
        $scriptsDir/conormaliseCanonicals.pl rm_$subject.fearful.wavered.1D rm_$subject.happy.wavered.1D rm_$subject.neutral.wavered.1D rm_$subject.sad.wavered.1D

        $scriptsDir/conormaliseCanonicals.pl \
    	    rm_$subject.fearful.remembered.wavered.1D rm_$subject.fearful.not.remembered.wavered.1D  rm_$subject.fearful.omission.wavered.1D \
    	    rm_$subject.happy.remembered.wavered.1D rm_$subject.happy.not.remembered.wavered.1D  rm_$subject.happy.omission.wavered.1D \
    	    rm_$subject.neutral.remembered.wavered.1D rm_$subject.neutral.not.remembered.wavered.1D  rm_$subject.neutral.omission.wavered.1D \
    	    rm_$subject.sad.remembered.wavered.1D rm_$subject.sad.not.remembered.wavered.1D  rm_$subject.sad.omission.wavered.1D

    	#echo "$subject: exists"
	
    	echo "Converting stim_times instruction regressors to binary 1D files"
    	for instr in A B C D ; do 
    	    timing_tool.py -timing ${subject}.acquisition.onsetDuration.instr${instr}.tab -timing_to_1D ${subject}.acquisition.onsetDuration.instr${instr}.binary.1D -tr 2 -min_frac 0.6 -run_len 860
    	done

	echo "####################################################################################################"
	echo "### Now making PPI regressors"
	$scriptsDir/makePpiContrastRegressors.pl --subject=${subject}

    else
	nonExistantSubjects="$nonExistantSubjects $subject"
	#echo "$subject: DOES NOT EXIST"
    fi
done

echo "Regressors for the the following subjects do not exist:"
echo "${nonExistantSubjects:-None}"

echo "Removing rm_ from the beginning of the normalized files"
for f in $( \ls -1 rm_*normalized.1D ) ; do
    mv -f $f $( echo $f | sed 's/rm_//g' ) 
done

find ./ -name "rm_*" -exec rm -f {} \;

#if [ ! -d ../Group.data/ ] ; then
#    mkdir ../Group.data/
#fi
##mv -f summary.tab ../Group.data/