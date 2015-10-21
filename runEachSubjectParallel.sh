#!/bin/bash

#set -x

##subjects=$( cd /Volumes/PROMISEPEGASUS/yangdata/cPine/data; \ls -1d *_[CD] )

####################################################################################################
# No resting state data for:
# 148_A 110_A 316_A 332_A 340_A
####################################################################################################
# No VBM anat data for:
# 341_A 346_A 347_A 348_A 354_A 357_A 365_A 367_A 343_A 358_A 359_A 360_A

#ctrlSubjects="107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 354_A 357_A 365_A 367_A"

#mddSubjects="110_A 111_A 112_A 113_A 114_A 117_A 134_A 144_A 147_A 158_A 160_A 161_A 300_A 301_A 304_A 313_A 316_A 323_A 324_A 336_A 340_A 343_A 359_A 137_A 310_A 322_A 325_A 344_A"

#mddSubjects="330_A 332_A 358_A"

##subjects="$ctrlSubjects $mddSubjects"

#subjects="168_A 350_A 369_A 374_A" 
##105_A

## controls only
#ctrl_subjects="107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 354_A 357_A 367_A 369_A 374_A 377_A" 

## MDD subjects only
#mdd_subjects="111_A 112_A 114_A 132_A 134_A 144_A 147_A 158_A 160_A 161_A 164_A 300_A 304_A 313_A 316_A 317_A 318_A 319_A 320_A 323_A 324_A 330_A 331_A 332_A 336_A 339_A 343_A 349_A 356_A 358_A 359_A 360_A 361_A 362_A 363_A 366_A"

#subjects="$ctrl_subjects $mdd_subjects"

#subjects="318_A 319_A 320_A 332_A"
##317_A 318_A 319_A 320_A 331_A 363_A 371_A 372_A 373_A 376_A"

subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.subjectList.txt )"


date

##for timepoint in A C D ; do

for timepoint in A  ; do

    for subject in $subjects ; do

	subject="${subject%_*}_${timepoint}"

	if [ -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subject/${subject}Pine+orig.HEAD ] ; then
	    
	# if [ ! -d $MDD_ROOT/log ] ; then 
	#     mkdir $MDD_ROOT/log 
	# fi
	    
	    cat <<EOF > run/run-$subject.sh
#!/bin/bash

AFNI_COMPRESSOR=GZIP
AFNI_DECONFLICT=OVERWRITE

if [ \$( uname -s ) == "Darwin" ] ; then
    CPUS=\$( sysctl -n hw.ncpu )
    MDD_ROOT=/Volumes/PROMISEPEGASUS/yangdata/cPine/
else 
## assume Linux
    CPUS=\$( expr \$( cat /proc/cpuinfo | grep -c processor ) - 4 )
    MDD_ROOT=/mnt/nfs/yangdata/cPine
fi

CPUS=4
JOBS="-jobs \$CPUS"
OMP_NUM_THREADS=\$CPUS

export CPUS JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT MDD_ROOT
#contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"
#contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad"
#contrasts="happyVsNeutral happyVsSad neutralVsSad"
#contrasts=fearfulVsHappyExtractedRightSgAcc

contrasts="happyVsNeutral happyVsSad neutralVsSad"


cd \$MDD_ROOT/scripts

# ./01-preprocessAnatomy.sh             --subject=$subject                               >  \$MDD_ROOT/log/$subject.log 2>&1
# ./01a-preprocessSubcorticalAnatomy.sh -s $subject                                      >  \$MDD_ROOT/log/$subject.first.log 2>&1
# ./02-preprocessFunctional.sh          --subject=$subject --drop=0 --fwhm=4.2 -c 0.3    >> \$MDD_ROOT/log/$subject.log 2>&1

#./02a-regenerateCensorFiles.sh --subject=$subject --drop=0 -c 0.3          > \$MDD_ROOT/log/$subject.log 2>&1

if [ ! -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt ] ; then 
#	./03-deconvolve-pine-withInstructionsRegressors.sh    --subject=$subject --polort=2       > \$MDD_ROOT/log/$subject.log 2>&1
#	./04-percentChange-pine-withInstructionsRegressors.sh --subject=$subject                 >> \$MDD_ROOT/log/$subject.log 2>&1

#	./03-deconvolve-pine.sh         --subject=$subject --polort=2       >> \$MDD_ROOT/log/$subject.log 2>&1
#	./03-deconvolve-pine-tent.sh    --subject=$subject --polort=2       >> \$MDD_ROOT/log/$subject.log 2>&1

#	./04-percentChange-pine.sh      --subject=$subject                  >> \$MDD_ROOT/log/$subject.log 2>&1

	for contrast in \$contrasts ; do
		./08-splitPpi.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/\${contrast}.seeds.txt -p 1 > \$MDD_ROOT/log/$subject.\$contrast.splitPpi.log 2>&1
	done

	# for mask in apriori.rois.withVmpfc.3mm+tlrc.HEAD ; do
	# 	maskName=\${mask%%+*}
	#  	for contrast in \$contrasts ; do
	# # 		# if [ -d /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subjectNumber/functional/ppi ] ; then 
	# # 		# 	mv -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subjectNumber/functional/ppi /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subjectNumber/functional/ppi-old
	# # 		# fi
	


	#  		./08-ppi.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/\${contrast}.\${maskName}.seeds.txt -p 1 -x \${maskName} > \$MDD_ROOT/log/$subject.\$contrast.\$maskName.ppi.log 2>&1
	
	# # 		./08-ppi.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/\$contrast.seeds.txt -p 1  > \$MDD_ROOT/log/$subject.\$contrast.ppi.log 2>&1
	# # 		./08-ppi.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt   -p 1 >> \$MDD_ROOT/log/$subject.\$contrast.apriori.ppi.log 2>&1
	
	# #		./extractPpiTimeseries.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/\$contrast.seeds.txt >> \$MDD_ROOT/log/$subject.\$contrast.ppi.log 2>&1
	# #		./extractPpiTimeseries.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt   >> \$MDD_ROOT/log/$subject.\$contrast.ppi.log 2>&1
	#  	done	for contrast in \$contrasts ; do
	#	./08-splitPpi.sh -s $subject  -c \$contrast -l ../data/config/ppi_seeds/\${contrast}.seeds.txt -p 1
	# done
	# done	

	#	./09-ppi-clean.sh --subject=$subject                  > \$MDD_ROOT/log/$subject.log 2>&1
	
	#	./01a-preprocessSubcorticalAnatomy.sh -s $subject > \$MDD_ROOT/log/$subject.first.log 2>&1
	#	./nonlinearfit.sh -s $subject > \$MDD_ROOT/log/$subject.log 2>&1
	#	./percentAUC.sh   -s $subject > \$MDD_ROOT/log/$subject.log 2>&1

else
	echo "*** WARNING: Found /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subject/functional/00_DO_NOT_ANALYSE_${subject}.txt. Skipping ${subject}."
fi

EOF
	    chmod +x  run/run-$subject.sh
	    echo $subject
	    nice -n 20 sem --no-notice --id pine -j1 ./run/run-$subject.sh ";" echo "$subject done"
	else
	    echo "*** /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subject/${subject}Pine+orig.HEAD does not exist. Skipping" 
	fi
    done
done
sem --no-notice --id pine --wait

date
