#!/bin/bash

#set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/data/sanDiego/cPine}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

GETOPT_OPTIONS=$( $GETOPT  -o "s:,f:" --longoptions "subject::" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided."
    exit
fi

cd $DATA/$subjectNumber 

if [ ! -f ${subjectNumber}.anat.nii.gz ] ; then
    if [ -f ${subjectNumber}_clp+orig.HEAD ] ; then 
	echo "*** NIfTI version of anatomy does not exist. Creating from clipped anatomy."
	3dcopy ${subjectNumber}_clp+orig.HEAD  ${subjectNumber}.anat.nii
    else
	echo "*** NIfTI version of anatomy does not exist. Creating."
	3dcopy ${subjectNumber}+orig.HEAD  ${subjectNumber}.anat.nii
    fi
fi

if [ ! -d anat ] ; then 
    mkdir anat
fi

if [ ! -f anat/${subjectNumber}.std.2.MNI.warpcoef.nii.gz ] ; then 

## anat_struc_brain.std.2.MNI.warpcoef.nii.gz ] ; then 
    
    echo "*** Reorienting NIfTI to standard"
    
    # if [ -f $subjectNumber.anat_clp.nii.gz ] ; then 
    #     echo "*** USING CLIPPED ANATOMY ***"
    # 	cp  $subjectNumber.anat_clp.nii.gz anat/$subjectNumber.anat.nii.gz
    # else
    # 	if [ ! -f $subjectNumber.anat.nii.gz ] ; then 
    # 	    3dcopy $subjectNumber+orig.HEAD $subjectNumber.anat.nii
    # 	fi
    # 	( cd anat; ln -sf ../$subjectNumber.anat.nii.gz . )
    # fi
    
    cd anat
    ln -sf ../$subjectNumber.anat.nii.gz .

    3dSkullStrip -input $subjectNumber.anat.nii.gz -prefix $subjectNumber.anat_struc_brain.nii
    3dcalc -a $subjectNumber.anat_struc_brain.nii.gz -expr 'step(a)' -prefix $subjectNumber.anat_struc_brain.mask.nii

    3dresample -orient RPI -inset $subjectNumber.anat_struc_brain.nii.gz -prefix $subjectNumber.anat_struc_brain.std.nii
    3dresample -orient RPI -inset $subjectNumber.anat.nii.gz -prefix $subjectNumber.anat_struc.std.nii
    3dresample -orient RPI -inset $subjectNumber.anat_struc_brain.mask.nii.gz -prefix $subjectNumber.anat_struc_brain.mask.std.nii

    echo "*** Computing affine alignment of skullstripped anat to MNI152 using FLIRT"
    flirt \
	-ref $MDD_STANDARD/MNI152_T1_3mm_brain.nii.gz \
	-in ${subjectNumber}.anat_struc_brain.std \
	-out ${subjectNumber}.anat_struc_brain.std.2.MNI.affine \
	-omat ${subjectNumber}.anat_struc_brain.std.2.MNI.affine.mat
    
    echo "*** Computing nonlinear alignment of anat to MNI152 using FNIRT"
    fnirt \
	--in=${subjectNumber}.anat_struc.std \
	--aff=${subjectNumber}.anat_struc_brain.std.2.MNI.affine.mat \
	--cout=${subjectNumber}.std.2.MNI.warpcoef \
	--iout=${subjectNumber}.std.2.MNI.nonlinear \
	--config=$MDD_STANDARD/T1_2_MNI152_3mm.cnf
    
    echo "*** Applying nonlinear warping skull stripped anatomy"
    applywarp \
	--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
	--in=${subjectNumber}.anat_struc_brain.std \
	--warp=${subjectNumber}.std.2.MNI.warpcoef \
	--out=${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear
    
    3dcopy ${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear.nii.gz ${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear
    
    echo "*** Generating mask from nonlinear warped skull stripped anatomy"
    3dAutomask -dilate 1 -prefix ${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear.mask ${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear.nii.gz

    echo "Segmenting brain for ${subjectNumber}"
    fast -t 1 -g -p -o segment ${subjectNumber}.anat_struc_brain.std.nii.gz 
    
    echo "Eroding WM mask by 1 Voxel"
    3dcalc -a segment_seg_2.nii.gz -prefix wm.eroded.nii \
	-b a+i -c a-i -d a+j -e a-j -f a+k -g a-k     \
	-expr 'a*(1-amongst(0,b,c,d,e,f,g))' 
    
    echo "Warping eroded WM mask to MNI space"
    applywarp \
	--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
	--in=wm.eroded \
	--warp=${subjectNumber}.std.2.MNI.warpcoef \
	--out=wm.eroded.MNI \
	--interp=nn
    
    echo "Finding overlap with WM prior"
    3dcalc -a $MDD_TISSUEPRIORS/3mm/avg152T1_white_bin.nii.gz -b wm.eroded.MNI.nii.gz -expr 'a*b' -prefix wm.eroded.masked.MNI
    3drefit -space MNI -view tlrc wm.eroded.masked.MNI+tlrc.HEAD
    
    echo "Warping eroded CSF mask to MNI space"
    applywarp \
	--ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
	--in=segment_seg_0.nii.gz \
	--warp=${subjectNumber}.std.2.MNI.warpcoef \
	--out=csf.MNI \
	--interp=nn
    
    echo "Finding overlap with CSF prior"
    3dcalc -a $MDD_TISSUEPRIORS/3mm/avg152T1_csf_bin.nii.gz -b csf.MNI.nii.gz -expr 'a*b' -prefix csf.masked.MNI
    3drefit -space MNI -view tlrc csf.masked.MNI+tlrc.HEAD
else
    echo "*** anat/${subjectNumber}.std.2.MNI.warpcoef.nii.gz already exists. Skipping this subject."
    echo "*** To force processing of ${subjectNumber} anatomy files, delete $DATA/$subjectNumber/anat/${subjectNumber}.std.2.MNI.warpcoef.nii.gz"

fi
