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

CPUS=2
JOBS="-jobs $CPUS"
OMP_NUM_THREADS=$CPUS

export CPUS JOBS OMP_NUM_THREADS AFNI_COMPRESSOR AFNI_DECONFLICT MDD_ROOT
contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad happyVsNeutral happyVsSad neutralVsSad"

#contrasts="fearfulVsHappy fearfulVsNeutral fearfulVsSad"

#contrasts="happyVsNeutral happyVsSad neutralVsSad"

cd $MDD_ROOT/scripts

# ./01-preprocessAnatomy.sh             --subject=132_A                               >  $MDD_ROOT/log/132_A.log 2>&1
# ./01a-preprocessSubcorticalAnatomy.sh -s 132_A                                      >  $MDD_ROOT/log/132_A.first.log 2>&1
# ./02-preprocessFunctional.sh          --subject=132_A --drop=0 --fwhm=4.2 -c 0.3    >> $MDD_ROOT/log/132_A.log 2>&1

#./02a-regenerateCensorFiles.sh --subject=132_A --drop=0 -c 0.3          > $MDD_ROOT/log/132_A.log 2>&1

if [ ! -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/132_A/functional/00_DO_NOT_ANALYSE_132_A.txt ] ; then 
#	./03-deconvolve-pine-withInstructionsRegressors.sh    --subject=132_A --polort=2       > $MDD_ROOT/log/132_A.log 2>&1
#	./04-percentChange-pine-withInstructionsRegressors.sh --subject=132_A                 >> $MDD_ROOT/log/132_A.log 2>&1

#	./03-deconvolve-pine.sh         --subject=132_A --polort=2       >> $MDD_ROOT/log/132_A.log 2>&1
#	./03-deconvolve-pine-tent.sh    --subject=132_A --polort=2       >> $MDD_ROOT/log/132_A.log 2>&1

#	./04-percentChange-pine.sh      --subject=132_A                  >> $MDD_ROOT/log/132_A.log 2>&1

#	./08a-ppi-clean.sh --subject=132_A                  >> $MDD_ROOT/log/132_A.log 2>&1


 	for contrast in $contrasts ; do
# 		# if [ -d /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi ] ; then 
# 		# 	mv -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi /Volumes/PROMISEPEGASUS/yangdata/cPine/data//functional/ppi-old
# 		# fi

# 		./08-ppi.sh -s 132_A  -c $contrast -l ../data/config/ppi_seeds/$contrast.seeds.txt -p 1 > $MDD_ROOT/log/132_A.$contrast.ppi.log 2>&1
# 		./08-ppi.sh -s 132_A  -c $contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt -p 1 > $MDD_ROOT/log/132_A.$contrast.apriori.ppi.log 2>&1

#		./extractPpiTimeseries.sh -s 132_A  -c $contrast -l ../data/config/ppi_seeds/$contrast.seeds.txt >> $MDD_ROOT/log/132_A.$contrast.ppi.log 2>&1
		./extractPpiTimeseries.sh -s 132_A  -c $contrast -l ../data/config/ppi_seeds/a-priori-seeds.txt   >> $MDD_ROOT/log/132_A.$contrast.ppi.log 2>&1
 	done

#	./01a-preprocessSubcorticalAnatomy.sh -s 132_A > $MDD_ROOT/log/132_A.first.log 2>&1

#	./nonlinearfit.sh -s 132_A > $MDD_ROOT/log/132_A.log 2>&1
#	./percentAUC.sh   -s 132_A > $MDD_ROOT/log/132_A.log 2>&1

else
	echo "*** WARNING: Found /Volumes/PROMISEPEGASUS/yangdata/cPine/data/132_A/functional/00_DO_NOT_ANALYSE_132_A.txt. Skipping 132_A."
fi

