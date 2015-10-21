#!/bin/bash

ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
DATA=$ROOT/data

subjects="$( cat ../data/config/control.subjectList.txt ../data/config/mdd.subjectList.txt )"

for subjectNumber in $subjects ; do
    cd  $DATA/$subjectNumber/functional/
    rm -fr ppi 
    rm -fr ppi-withCleanUp/ppi
    
    #rm -f *.contrast.1D *.contrast.detrended.1D  *.contrast.detrended.transposed.1D  hrf.1D  *.contrast.detrended.transposed.convolved.1D  *.err
    #rm -f *.contrast.interactionTerm.1D *.contrast.interactionTerm.wavered.1D  *.dec+tlrc.*  *xmat* *REML_cmd  *z-score* *contrast.detrended.transposed.convolved.transposed.1D

done