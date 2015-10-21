#!/bin/bash

MDD_ROOT=/Volumes/PROMISEPEGASUS/yangdata/cPine/

cd ../data/Group.data

#echo "subject,seed,contrast"
for ff in nodata*mpfc* ; do 
    badSubjectCount=$( cat $ff | wc -l )
    
    badSubjects=( $( cat $ff ) ) 
    
    #seedInfo="$( echo $ff | gawk -F '.' '{print $3 "," $4 "." $5 "." $6 "." $7;}')"
    seedInfo="$( echo $ff | gawk -F '.' '{print $3 "," $5 "." $6 "." $7 "." $8;}')"
    
    roi="$( echo $ff | gawk -F '.' '{print $3 ;}' | sed "s/roi//" )"
    contrast="$( echo $ff | gawk -F '.' '{print $5 }' | sed "s/apriori//" )"
    contrastAndMask="$( echo $ff | gawk -F '.' '{print $5 "." $6 "." $7 "." $8;}')"
    mask=$( echo $contrastAndMask | sed "s/$contrast//" )

    for (( ii=0; ii < $badSubjectCount; ii=ii+1 ))
    do
	subject=${badSubjects[$ii]}
	if [ $subject != 137_A ] && [ $subject != 121_A ] ; then
	    echo "Running script for subject: ${subject} roi: $roi contrast: $contrast mask: $mask"
	    ( cd $MDD_ROOT/scripts; ./08-ppi.sh -s ${subject} -c $contrast -l ../data/config/ppi_seeds/${contrast}.${mask}.seeds.txt -p 1 -x ${mask} -b $roi > $MDD_ROOT/log/$subject.$contrast.$mask.roi$roi.ppi.log 2>&1 )
	fi

    done
done
