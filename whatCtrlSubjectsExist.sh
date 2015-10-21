#!/bin/bash

#set -x 

dataRoot="/Volumes/PROMISEPEGASUS/yangdata/"

## my original list 
## 107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 354_A 357_A 365_A 367_A

## new list
newCtrlList="105_A 107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 156_A 157_A 159_A 162_A 163_A 165_A 168_A 302_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 334_A 337_A 341_A 346_A 347_A 348_A 350_A 354_A 357_A 365_A 367_A 369_A 374_A 377_A 379_A 380_A"

for subject in $newCtrlList ; do

    line="$subject"
    if [ -d $dataRoot/$subject ] ; then
	line="$line YES"
    else 
	line="$line NO"
    fi

    if [ -d $dataRoot/cPine/data/$subject ] ; then
	line="$line YES"
    else 
	line="$line NO"
    fi

    if [ -f $dataRoot/cPine/data/$subject/${subject}Pine+orig.HEAD ] ; then
	line="$line YES"
    else 
	line="$line NO"
    fi
    echo $line
done
