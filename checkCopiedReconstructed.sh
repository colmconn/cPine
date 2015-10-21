#!/bin/bash

mdd_subjects="111_A 112_A 114_A 132_A 134_A 144_A 147_A 158_A 160_A 161_A 164_A 300_A       304_A 313_A 316_A 317_A 318_A 319_A 320_A 323_A 324_A 330_A 331_A 332_A                   336_A 339_A 340_A 343_A 349_A 356_A 358_A 359_A 360_A 361_A 362_A 363_A 366_A 371_A 372_A 373_A 376_A"

ctrl_subjects="107_A 108_A 109_A 116_A 119_A 121_A 122_A 123_A 124_A 126_A 127_A 131_A 135_A 138_A 139_A 141_A 142_A 143_A 145_A 146_A 148_A 151_A 152_A 153_A 155_A 157_A 159_A 162_A 163_A 165_A 303_A 306_A 307_A 308_A 312_A 321_A 326_A 328_A 337_A 341_A 346_A 347_A 348_A 354_A 357_A 367_A 369_A 374_A 375_A 377_A 380_A 382_A 386_A"


cd /mnt/nfs/yangdata/

header="subject,group,reconNotInCPine,reconInCPine,notRecon"

function checkExistance {

    local group=$1
    local subjects=$2

    for f in $subjects ; do 

	line="$f,$group"
	if [ -f /Volumes/PROMISEPEGASUS/yangdata/${f}BRIKS/${f}Pine+orig.HEAD ] && [ ! -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/${f}/${f}Pine+orig.HEAD ] ; then
	    line="${line},TRUE"
	else 
	    line="${line},FALSE"
	fi

	if [ -f /Volumes/PROMISEPEGASUS/yangdata/${f}BRIKS/${f}Pine+orig.HEAD ] && [ -f /Volumes/PROMISEPEGASUS/yangdata/cPine/data/${f}/${f}Pine+orig.HEAD ] ; then
	    line="${line},TRUE"
	else 
	    line="${line},FALSE"
	fi

	if [ ! -f /Volumes/PROMISEPEGASUS/yangdata/${f}BRIKS/${f}Pine+orig.HEAD ] ; then
	    line="${line},TRUE"
	else 
	    line="${line},FALSE"
	fi
	echo $line

    done

}

echo $header
checkExistance mdd "$mdd_subjects"
checkExistance ctrl "$ctrl_subjects"




