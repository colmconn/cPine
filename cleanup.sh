#!/bin/bash

#set -x

subjects=$( cd /Volumes/PROMISEPEGASUS/yangdata/cPine/data; \ls -1d *_[ACD] )

for subject in $subjects ; do

    echo "Cleaning $subject"
    cd /Volumes/PROMISEPEGASUS/yangdata/cPine/data/$subject/functional && rm -f  *{dec,fitts,errts,%cs,xmat,jpg,err,contrast,REML}*

done