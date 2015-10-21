#!/bin/bash

#set -x 

##subjects="389_A 391_A 392_A 393_A 395_A 396_A 398_A 399_A 400_A 401_A 402_A 403_A 404_A 505_A"

#subjects="379_A 380_A 381_A 382_A 383_A 384_A 385_A 386_A 387_A 388_A 389_A 390_A 392_A 393_A 394_A 397_A 405_A 406_A 407_A"

subjects="386_A"

for subject in $subjects ; do 

    for f in ESTOP_Results.txt  FRT_Results.txt Reap_resultsmri.txt ResultsMorhper.txt
    do 
	subject=${subject%%_A}
	if grep -q "^10$subject" /Volumes/PROMISEPEGASUS/yangdata/404_A/404_A_beh/$f
	then
	    echo "Found $subject in $f"
	else 
	    echo "$subject is *NOT* in $f"
	fi
    done
done