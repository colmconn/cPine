#!/bin/bash
stimuli="fearful happy"

echo $stimuli
arr=($(echo $stimuli | tr " " "\n" ))

#$( echo ${arr[1]}[0:1] | tr '[:lower;]' '[:upper]' )${arr[1]}[1:
cr="${arr[0]}Vs$(echo ${arr[1]} | gsed 's/.*/\u&/')"

echo $cr

