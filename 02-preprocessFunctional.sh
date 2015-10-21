#!/bin/bash

set -x 

# if ctrl-c is typed exit immediatly
trap exit SIGHUP SIGINT SIGTERM

programName=`basename $0`

GETOPT=$( which getopt )
ROOT=${MDD_ROOT:-/Volumes/PROMISEPEGASUS/yangdata/cPine}
DATA=$ROOT/data
MDD_STANDARD=$ROOT/standard
MDD_TISSUEPRIORS=$ROOT/tissuepriors
scriptsDir=${ROOT}/scripts

task=pine
censorThreshold=0.2

GETOPT_OPTIONS=$( $GETOPT  -o "s:d:f:c:" --longoptions "subject:,drop:,fwhm:,censorThreshold:" -n ${programName} -- "$@" )
exitStatus=$?
if [ $exitStatus != 0 ] ; then 
    echo "Error with getopt. Terminating..." >&2 
    exit $exitStatus
fi

drop=0

# Note the quotes around `$GETOPT_OPTIONS': they are essential!
eval set -- "$GETOPT_OPTIONS"
while true ; do 
    case "$1" in
	-s|--subject)
	    subjectNumber=$2; shift 2 ;;
	-d|--drop)
	    drop=$2; shift 2 ;;
#	-p|--polort)
#	    polort=$2; shift 2 ;;
	-f|--fwhm)
	    fwhm=$2 ;
	    sigma=$( echo "scale=10 ; ${fwhm}/2.3548" | bc ) ;
	    shift 2 ;;
	-c|--censorThreshold)
	    censorThreshold=$2; shift 2 ;;
	--) 
	    shift ; break ;;

	*) 
	    echo "${programName}: ${1}: invalid option" >&2
	    exit 2 ;;
    esac
done

echo "Using $censorThreshold as threshold of censoring motion"

if [ -z $subjectNumber ] ; then 
    echo "*** ERROR: The subject ID was not provided."
    exit
fi

if [ -z $fwhm ] ; then 
    echo "*** ERROR: The fwhm was not provided."
    exit
fi

prefixPrefix=${subjectNumber}.pine

cd  $DATA/$subjectNumber/
mkdir functional
cd functional

##if [ 1 == 0 ] ; then

if [ ! -f ${subjectNumber}Pine+orig.HEAD ] ; then 
    cp ../${subjectNumber}Pine+orig.HEAD  .
    cp ../${subjectNumber}Pine+orig.BRIK* .
fi

## set the slice time offsets in the attributes as to3d can't extract
## them from the DICOMS properly (or they are not in the DICOMS in the
## first instance)

tr=$( 3dinfo -tr ${subjectNumber}Pine+orig.HEAD ) 
if [ $tr == "2000" ] ; then 
    3drefit -atrfloat TAXIS_OFFSETS '0 1000 50 1050 100 1100 150 1150 200 1200 250 1250 300 1300 350 1350 400 1400 450 1450 500 1500 550 1550 600 1600 650 1650 700 1700 750 1750 800 1800 850 1850 900 1900 950 1950' ${subjectNumber}Pine+orig.HEAD
else 
    3drefit -atrfloat TAXIS_OFFSETS '0 1 0.05 1.05 0.1 1.1 0.15 1.15 0.2 1.2 0.25 1.25 0.3 1.3 0.35 1.35 0.4 1.4 0.45 1.45 0.5 1.5 0.55 1.55 0.6 1.6 0.65 1.65 0.7 1.7 0.75 1.75 0.8 1.8 0.85 1.85 0.9 1.9 0.95 1.95' ${subjectNumber}Pine+orig.HEAD
fi

echo "*** Generating outlier count"
3dToutcount -automask -fraction -polort 4 -legendre \
    ${subjectNumber}Pine+orig.HEAD > outcount.1D

# censor outlier TRs per run, ignoring the first 0 TRs
# - censor when more than 0.1 of automask voxels are outliers
# - step() defines which TRs to remove via censoring
1deval -a outcount.1D -expr "1-step(a-0.1)" > rm.outliers.censor.1D

#echo "*** Deoblquing resting state data"
#3dWarp -deoblique -prefix ${prefixPrefix} ${subjectNumber}Pine+orig.HEAD
3dcopy ${subjectNumber}Pine+orig.HEAD ${prefixPrefix}
rm -f ${subjectNumber}Pine+orig.*

cp ../anat/$subjectNumber.anat_struc_brain.std.nii.gz  .
cp ../anat/$subjectNumber.anat_struc_brain.nii.gz  .
3dcopy $subjectNumber.anat_struc_brain.std.nii.gz $subjectNumber.anat_struc_brain.std
3dcopy $subjectNumber.anat_struc_brain.nii.gz $subjectNumber.anat_struc_brain

echo "*** Aligning"
align_epi_anat.py -anat $subjectNumber.anat_struc_brain+orig.HEAD \
    -tshift on -tshift_opts "-Fourier" \
    -anat_has_skull no -epi_strip 3dAutomask \
    -epi ${prefixPrefix}+orig.HEAD -epi_base 0 \
    -epi2anat -ex_mode echo

echo "*** Reorienting to RPI"
3dresample -orient RPI -inset ${prefixPrefix}_al+orig.HEAD -prefix ${prefixPrefix}.std_al

if [ ! -f ${prefixPrefix}.std_al+orig.HEAD ] ; then
    echo "***ERROR: ${prefixPrefix}.std_al+orig.HEAD does not exist. Something went wrong with the alignment. Exiting."
    exit
fi

if [ -f $subjectNumber.pine_tsh_vr_motion.1D ] ; then 
    motionFile=$subjectNumber.pine_tsh_vr_motion.1D
else
    motionFile=$subjectNumber.pine_vr_motion.1D
fi

if [[ $drop != 0 ]] ; then 
    echo "*** Dropping first TRs"
    3dcalc -a ${prefixPrefix}.std_al+orig.HEAD"[$drop..$]" -expr 'a' -prefix ${prefixPrefix}.trimmed
    tail -n +$( expr $drop + 1 ) $motionFile > ${prefixPrefix}.motion.1D
    tail -n +$( expr $drop + 1 ) rm.outliers.censor.1D > rm.outliers.censor.trimmed.1D
else
    echo "*** Not dropping any volumes as drop = $drop"
    cat $motionFile > ${prefixPrefix}.motion.1D
fi

# compute de-meaned motion parameters (for use in regression)
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -demean -write motion_demean.1D

# compute motion parameter derivatives (for use in regression)
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -derivative -demean -write motion_deriv.1D

# create censor file motion_${subj}_censor.1D, for censoring motion 
1d_tool.py -overwrite -infile ${prefixPrefix}.motion.1D -set_nruns 1                              \
    -set_tr $tr -show_censor_count -censor_prev_TR                           \
    -censor_motion ${censorThreshold} motion_${subjectNumber}

if [[ $drop != 0 ]] ; then
## now combine both the motion and outlier censor files into one
    1deval -a motion_${subjectNumber}_censor.1D -b rm.outliers.censor.trimmed.1D  \
	-expr "a*b" > censor_${subjectNumber}_combined_2.1D
else 
    1deval -a motion_${subjectNumber}_censor.1D -b rm.outliers.censor.1D  \
	-expr "a*b" > censor_${subjectNumber}_combined_2.1D
fi

echo "*** Smoothing and masking"
#fslmaths ${prefixPrefix}.masked.nii.gz -kernel gauss ${sigma} -fmean -mas ${prefixPrefix}.mask.nii.gz ${prefixPrefix}.smoothed

if [[ $drop != 0 ]] ; then
    3dresample -master ${prefixPrefix}.trimmed+orig.HEAD \
	-prefix ./${subjectNumber}.anat_struc_brain.mask.std.resampled \
	-inset ../anat/${subjectNumber}.anat_struc_brain.mask.std.nii.gz
    3dBlurInMask -mask ./${subjectNumber}.anat_struc_brain.mask.std.resampled+orig.HEAD -FWHM $fwhm -prefix ${prefixPrefix}.smoothed ${prefixPrefix}.trimmed+orig.HEAD

else 
    3dresample -master ${prefixPrefix}.std_al+orig.HEAD \
	-prefix ./${subjectNumber}.anat_struc_brain.mask.std.resampled \
	-inset ../anat/${subjectNumber}.anat_struc_brain.mask.std.nii.gz
    3dBlurInMask -mask ./${subjectNumber}.anat_struc_brain.mask.std.resampled+orig.HEAD -FWHM $fwhm -prefix ${prefixPrefix}.smoothed ${prefixPrefix}.std_al+orig.HEAD
fi

##3dBlurInMask -mask ./${subjectNumber}.anat_struc_brain.mask.std.resampled+orig.HEAD -FWHM $fwhm -prefix ${prefixPrefix}.smoothed ${prefixPrefix}.trimmed+orig.HEAD
#3dBlurInMask -mask ./${subjectNumber}.anat_struc_brain.mask.std.resampled+orig.HEAD -FWHM $fwhm -prefix ${prefixPrefix}.preprocessed.nii ${prefixPrefix}.trimmed+orig.HEAD

echo "*** Grand mean scaling"
3dTstat -mean -prefix ${prefixPrefix}.mean ${prefixPrefix}.smoothed+orig.HEAD
3dcalc -float -a ${prefixPrefix}.smoothed+orig.HEAD -b ${prefixPrefix}.mean+orig.HEAD  -expr "(a/b)*10000" -prefix ${prefixPrefix}.preprocessed

3dcopy ${prefixPrefix}.preprocessed+orig.HEAD ${prefixPrefix}.preprocessed.nii

echo "*** Creating mask from preprocessed data"
fslmaths ${prefixPrefix}.preprocessed.nii.gz -Tmin -bin ${prefixPrefix}.preprocessed.mask.nii.gz -odt char

3dbucket -prefix ${prefixPrefix}.func.native.std.nii ${prefixPrefix}.std_al+orig.HEAD\[0\] 

echo "*** Resampling ${prefixPrefix}.func.native.std.nii.gz to original space anatomy gridset"
3dresample -master ../anat/${subjectNumber}.anat_struc.std.nii.gz -prefix ./$prefixPrefix.func.std.resampled.nii \
    -inset ${prefixPrefix}.func.native.std.nii.gz

echo "*** Aligning $subjectNumber.${task}.func.native.std to anatomy"
flirt -in ${prefixPrefix}.func.native.std -ref ${prefixPrefix}.func.std.resampled \
    -out ${prefixPrefix}.func2anat.flirt -omat ${prefixPrefix}.func2anat.flirt.mat

echo "*** Applying nonlinear warp to MNI space to ${prefixPrefix}.preprocessed"
applywarp \
    --ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
    --in=${prefixPrefix}.preprocessed  \
    --warp=../anat/${subjectNumber}.std.2.MNI.warpcoef \
    --premat=${prefixPrefix}.func2anat.flirt.mat \
    --out=${prefixPrefix}.preprocessed.MNI

echo "*** Applying nonlinear warp to MNI space to ${prefixPrefix}.preprocessed.mask"
applywarp \
    --ref=$MDD_STANDARD/MNI152_T1_3mm.nii.gz \
    --in=${prefixPrefix}.preprocessed.mask \
    --warp=../anat/${subjectNumber}.std.2.MNI.warpcoef \
    --premat=${prefixPrefix}.func2anat.flirt.mat \
    --out=${prefixPrefix}.preprocessed.mask.MNI

##fi

3dcopy ${prefixPrefix}.preprocessed.mask.MNI.nii.gz ${prefixPrefix}.preprocessed.mask.MNI
3dcopy ${prefixPrefix}.preprocessed.MNI.nii.gz ${prefixPrefix}.preprocessed.MNI
 
## for some anoying the TR is not preserved by the applywarp so make
## sure it's set correctly here by using 3drefit
3drefit -space MNI -TR $tr ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD
3drefit -space MNI ${prefixPrefix}.preprocessed.mask.MNI+tlrc.HEAD

3dresample -master ${prefixPrefix}.preprocessed.MNI+tlrc.HEAD -inset ../anat/csf.masked.MNI+tlrc.HEAD       -prefix ./csf.masked.func.MNI
3dresample -master ${prefixPrefix}.preprocessed.MNI.nii.gz    -inset ../anat/wm.eroded.masked.MNI+tlrc.HEAD -prefix ./wm.eroded.masked.func.MNI

3drefit -space MNI -view tlrc ./csf.masked.func.MNI+tlrc.HEAD
3drefit -space MNI -view tlrc ./csf.masked.func.MNI+tlrc.HEAD

cp ../anat/${subjectNumber}.anat_struc_brain.std.2.MNI.nonlinear+tlrc.* ./

motionThresholdPrecentage=0.2
numberOfCensoredVolumes=$( cat censor_${subjectNumber}_combined_2.1D | gawk '{a+=(1-$0)}END{print a}' )
totalNumberOfVolumes=$( cat censor_${subjectNumber}_combined_2.1D | wc -l )
twentyPercent=$( echo "scale=0; $motionThresholdPrecentage*$totalNumberOfVolumes" | bc | cut -f 1 -d '.' )

if [[ $numberOfCensoredVolumes -gt $twentyPercent ]] ; then 
    echo "A total of $numberOfCensoredVolumes of $totalNumberOfVolumes we censored which is greater than 20% (n=$twentyPercent) of all total volumes of this subject" >  00_DO_NOT_ANALYSE_${subjectNumber}.txt
    echo "*** WARNING: This subject will not be analysed due to having more that 20% of their volumes censored."
fi