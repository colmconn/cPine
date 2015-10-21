#!/bin/bash

trap exit SIGHUP SIGINT SIGTERM

set -x 

subjects="$*"

for subject in $subjects ; do 
    if [ -d /Volumes/PROMISEPEGASUS/yangdata/$subject ] ; then 
	if [ ! -d /Volumes/PROMISEPEGASUS/yangdata/$subject/${subject}_BEH ] || [ ! -d /Volumes/PROMISEPEGASUS/yangdata/$subject/${subject}_beh ] ; then 
	    scp -r taisen:/mnt/nfs/yangdata/MDD_3TW/data/${subject}_BEH /Volumes/PROMISEPEGASUS/yangdata/$subject/
	    scp -r taisen:/mnt/nfs/yangdata/MDD_3TW/data/${subject}_beh /Volumes/PROMISEPEGASUS/yangdata/$subject/
	fi
    else 
	echo "*** No such directory:  /Volumes/PROMISEPEGASUS/yangdata/$subject"
	echo "Skipping"
    fi
done
