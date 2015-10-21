#!/bin/bash

AFNI_COMPRESSOR=GZIP
AFNI_DECONFLICT=OVERWRITE

if [ $( uname -s ) == "Darwin" ] ; then
    CPUS=$( sysctl -n hw.ncpu )
    MDD_ROOT=/Volumes/PROMISEPEGASUS/yangdata/cPine/
else 
## assume Linux
    CPUS=$( expr $( cat /proc/cpuinfo | grep -c processor ) - 4 )
    MDD_ROOT=/mnt/nfs/yangdata/cPine
fi

CPUS=4
JOBS="-jobs $CPUS"
OMP_NUM_THREADS=$CPUS

export CPUS JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT MDD_ROOT
#contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"
#contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad"
#contrasts="happyVsNeutral happyVsSad neutralVsSad"
#contrasts=fearfulVsHappyExtractedRightSgAcc

contrasts="happyVsNeutral happyVsSad neutralVsSad"


cd $MDD_ROOT/scripts

# ./01-preprocessAnatomy.sh             --subject=155_A                               >  $MDD_ROOT/log/155_A.log 2>&1
# ./01a-preprocessSubcorticalAnatomy.sh -s 155_A                                      >  $MDD_ROOT/log/155_A.first.log 2>&1
# ./02-preprocessFunctional.sh          --subject=155_A --drop=0 --fwhm=4.2 -c 0.3    >> $MDD_ROOT/log/155_A.log 2>&1

#./02a-regenerateCensorFiles.sh --subject=155_A --drop=0 -c 0.3          > $MDD_ROOT/log/155_A.log 2>&1

if [ ! -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/155_A/functional/00_DO_NOT_ANALYSE_155_A.txt ] ; then 
#	./03-deconvolve-pine-withInstructionsRegressors.sh    --subject=155_A --polort=2       > $MDD_ROOT/log/155_A.log 2>&1
#	./04-percentChange-pine-withInstructionsRegressors.sh --subject=155_A                 >> $MDD_ROOT/log/155_A.log 2>&1

#	./03-deconvolve-pine.sh         --subject=155_A --polort=2       >> $MDD_ROOT/log/155_A.log 2>&1
#	./03-deconvolve-pine-tent.sh    --subject=155_A --polort=2       >> $MDD_ROOT/log/155_A.log 2>&1

#	./04-percentChange-pine.sh      --subject=155_A                  >> $MDD_ROOT/log/155_A.log 2>&1

	for contrast in $contrasts ; do
		./08-splitPpi.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/${contrast}.seeds.txt -p 1 > $MDD_ROOT/log/155_A.$contrast.splitPpi.log 2>&1
	done

	# for mask in apriori.rois.withVmpfc.3mm+tlrc.HEAD ; do
	# 	maskName=${mask%%+*}
	#  	for contrast in $contrasts ; do
	# # 		# if [ -d /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi ] ; then 
	# # 		# 	mv -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi-old
	# # 		# fi
	


	#  		./08-ppi.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/${contrast}.${maskName}.seeds.txt -p 1 -x ${maskName} > $MDD_ROOT/log/155_A.$contrast.$maskName.ppi.log 2>&1
	
	# # 		./08-ppi.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/$contrast.seeds.txt -p 1  > $MDD_ROOT/log/155_A.$contrast.ppi.log 2>&1
	# # 		./08-ppi.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt   -p 1 >> $MDD_ROOT/log/155_A.$contrast.apriori.ppi.log 2>&1
	
	# #		./extractPpiTimeseries.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/$contrast.seeds.txt >> $MDD_ROOT/log/155_A.$contrast.ppi.log 2>&1
	# #		./extractPpiTimeseries.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt   >> $MDD_ROOT/log/155_A.$contrast.ppi.log 2>&1
	#  	done	for contrast in $contrasts ; do
	#	./08-splitPpi.sh -s 155_A  -c $contrast -l ../data/config/ppi_seeds/${contrast}.seeds.txt -p 1
	# done
	# done	

	#	./09-ppi-clean.sh --subject=155_A                  > $MDD_ROOT/log/155_A.log 2>&1
	
	#	./01a-preprocessSubcorticalAnatomy.sh -s 155_A > $MDD_ROOT/log/155_A.first.log 2>&1
	#	./nonlinearfit.sh -s 155_A > $MDD_ROOT/log/155_A.log 2>&1
	#	./percentAUC.sh   -s 155_A > $MDD_ROOT/log/155_A.log 2>&1

else
	echo "*** WARNING: Found /Volumes/PROMISEPEGASUS/yangdata/cPine/data/155_A/functional/00_DO_NOT_ANALYSE_155_A.txt. Skipping 155_A."
fi

