#!/bin/csh
#Set up FSL environmental variables
if ($#argv <12) then
echo "***FCProcess: A script for preprocessing and individual-level analysis of functional connectivity and psychophysiological interactions"
echo "--------------------------------------------------------------"
echo "With modifications according to Jo et al., 2010 (artifact detection and removal, local white-matter regressors, etc.)."
echo "--------------------------------------------------------------"
echo "If you find this script useful, please cite: Jo, H.J., Saad, Z.S., Simmons, W. K., Milbury, L. A., Cox, R.W.  Mapping sources of correlation in resting state FMRI, with artifact detection and removal.  Neuroimage, 2010, 52:2, 571-582."
echo " -------------------------------------------------------------"
echo "HELP NEEDED? READ THE INFO BELOW FIRST, THEN CONTACT GREG FONZO WITH REMAINING QUESTIONS AT GFONZO@UCSD.EDU OR 858-246-0622"
echo "--------------------------------------------------------------"
echo "Below is a list of helpful tips for using this script."
echo "--------------------------------------------------------------"
echo "The proper usage is: ./FCProcess [subj#] [taskname] [studydirectory] [scripttype] [skullstrip suffix] [talairach suffix] [talairach type] [motion regressor file] [blurring radius] [ROI mask list] [FCtype] [censor file name]"
echo "--------------------------------------------------------------"
echo "Definitions of terms:"
echo "############"
echo "subj# = The numerical identifier for the subject (e.g., 0716198301)"
echo "############"
echo "taskname = The name of the task you would like to analyze (e.g., hariri, stimex, etc.)  This should include any suffixes or prefixes which were added during the preprocessing steps (e.g., hariri_dtsali2anat), or basically everything after the subj# but before the suffix (+orig,+tlrc)."
echo "############"
echo "studydirectory = The study directory where your data is located, with full path info (e.g., /mnt/nfs/studies/INSULAfmri)"
echo "############"
echo "scripttype = There are two modes for this script.  Enter 1 if you would like to use it in overwrite mode.  This will erase any prior files which might have been created in previous runs of the script.  This is useful if you made some mistakes earlier on and would like to start over.  Enter 2 if you would like to keep files which have been already created in earlier runs of the script.  This is useful if you needed to quit out of the script for some reason but don't want to have to reprocess the files which have already been created."
echo "############"
echo "skullstrip suffix = The suffix placed after the subj# to denote a skull-stripped anatomical (e.g., SS, NS, noskull, etc.)"
echo "############"
echo "talairach suffix = The suffix placed after the subj# to denote a talairached anatomical (e.g., AutoTal, Tal, etc.).  This is different than the normal +tlrc suffix.  If you used the auto-talairaching script, you might have added a suffix to the talairached file after the subj# but before the +tlrc suffix.  If so, enter that here.  If you did not enter a suffix for the talairached anatomical (e.g., 0716198391+tlrc), then please create one for the purposes of this script by using the 3dcopy command (e.g., 3dcopy 0716198301+tlrc 0716198301Tal+tlrc)."
echo "############"
echo "talairach type = Enter 1 here if you used the auto-talairaching script (@auto_tlrc).  Enter 2 here if you used the time-consuming, manual, gui-based talairaching, and next time save yourself some time and energy and just use the auto-talairaching feature.  Life is too short to waste it trying to identify the anterior and posterior commissure.  Unless of course your research advisor tells you to do it the hard way, in which case you should probably listen."
echo "############"
echo "motion regressor file = The name of the file containing the motion regressors from the realignment stage, usually output from 3dvolreg.  I believe it follows the -dfile option on the command line.  It should be located in the subject's BRIKS directory."
echo "############"
echo "blurring radius = The numerical value (in millimeters) of the blurring radius you would like to use for the FC analysis (e.g., 4.0, 6.0)."
echo "############"
echo "ROI mask list = This is a text file with the names of the brik files which contain the seed regions for the FC analysis (e.g., MyROIforFC+tlrc).  The masks themselves should be placed in the testdir directory within your study directory, and the text file should be placed within your scripts directory.  Please only have 1 SEED REGION or SEED CLUSTER per brik file.  It should have a uniform numerical value (e.g., all voxels in cluster/seed set to 1) as output, for example, by the SaveMask command in the clustering window or using the -1clust_order option from 3dmerge.  Within the text file, on each line have the first argument be the name of the mask (e.g., MYROIforFC+tlrc), the second argument the name of the contrast your would like to model (e.g., Fear-Happy, Negative-Positive, no spaces between words), and the third argument the name of the seed region or cluster (e.g., LAntIns, RAmygdala, no spaces between words).  Your file can have multiple lines to process multiple seed regions/clusters."
echo "############"
echo "FCtype = Enter 1 here if you want to look at a psychophysiological interaction (PPI), i.e., context-dependent connectivity.  Enter 2 here if you want to look at functional connectivity across the entire timeseries without regard to task conditions.  IF YOU WOULD LIKE TO DO OPTION 2, YOU WILL NEED TO MODIFY THE 3DDECONVOLVE COMMAND AT THE END OF THE SCRIPT TO ADD IN THE NUISANCE REGRESSORS (E.G., TASK CONDITION REGRESSORS).  If you don't know how to do this or don't understand the difference between a PPI and a regular FC analysis, please seek help immediately!  I would recommend asking someone with lots of time on their hands (i.e., probably not the author of this script, unless you are willing to bribe him with sugary foods or cash money, in which case he will be happy to help you)."
echo "IF TYPE 1 (PPI), you will also need a one-column text file with 1's, 0's, and -1's corresponding to the timepoints of the contrast you would like to model (e.g., if I want to look at connectivity for matching fearful vs. happy faces, I would code all the timepoints (prior to convolution with the HRF) with fear faces 1, all the timepoints with happy faces -1, and all the rest 0).  Please name this file [taskname]FC[contrastname].1D and place it in your scripts directory.  You will also need a regressor file corresponding to the contrast condition convolved with the HRF to model in the 3dDeconvolve analysis.  Please name this [contrastname].1D and place it in your scripts directory."
echo "#############"
echo "censor file name = The name of the one-column file denoting the outlier timepoints you would like to remove from the deconvolution analysis.  The first 2-3 timepoints are almost always censored out due to the scanner reaching thermal equilibirum, as well as others relating to movement artifacts.  Please place this file in the subject's BRIKS directory."
echo "#############"
echo "If you move this file to your own directory, be sure to copy the following scripts from the original directory (/mnf/nfs/scripts): 3ddExtallrFz and 3ddExtallCoef. Otherwise, you'll receive errors when you try to run this script"
echo "############"
echo "That's all he wrote...Use your best judgment.  Read through the script first before using it so you have some idea what's going on in the background.  Or ask someone with more experience to help you.  If you ask nicely, they will usually oblige...usually.  Nothing is certain in this life.  But if you don't know that you shouldn't be in this line of research anyway.  Get out while you can!"
echo "IF YOU ARE SEEING THIS AND WERE EXPECTING TO RUN THE SCRIPT, YOU MESSED UP SOMEWHERE.  GO BACK AND CHECK YOUR WORK.  GOOD LUCK!"
	exit 1
endif
echo "FCProcess: A script for preprocessing and individual-level analysis of functional connectivity and psychophysiological interactions"
echo "With modifications according to Jo et al., 2010 (artifact detection and removal, local white-matter regressors, etc.)."
echo "If you find this script useful, please cite: Jo, H.J., Saad, Z.S., Simmons, W. K., Milbury, L. A., Cox, R.W.  Mapping sources of correlation in resting state FMRI, with artifact detection and removal.  Neuroimage, 2010, 52:2, 571-582."
echo "If you would like to see the help screen, please type ./FCProcess in the terminal window."
set subject = $1
set task = $2
#This should be the name of the preprocessed timeseries (e.g., motion-corrected, realigned, despiked, etc.)  It should include any suffixes or prefixes following the subject number and before the +orig (e.g., hariri_dtsali2anat).
set studydirectory = $3
#The study directory name with full path, e.g., /mnt/nfs/studies/INSULAfmri
set scripttype = $4
#Enter '1' for overwrite mode, enter '2' not to erase files which have been previously written
set SSsuffix = $5
#This should be the suffix after the subject # which is used to denote a skull-stripped anatomical, e.g., SS
set talsuffix = $6
#This should be the suffix after the subject # which is used to denote a talairached dataset.  This is useful if you used the auto-taliarching script where sometimes people add a suffix to the dataset name, e.g., SS1.0AutoTal.  If you don't have a suffix after your talairached dataset (e.g., 1111201201+tlrc), then please give it one for the purposes of this script.  You can use the 3dcopy command to do this (e.g., 3dcopy 1111201201+tlrc 1111201201Tal+tlrc).
set taltype = $7
#Enter '1' if you used the auto-talairaching program, enter '2' if you used manual talairaching
set regfile = $8
#This is the name of the movement regressor file output from 3dvolreg, usually denoted in the command line following the -d option
set blur = $9
#Enter a numerical value, e.g., 4.0 or 6.0
set roimasklist = $10
#This should be a text file with 3 columns and multiple lines, if needed; 1st column is the name of the ROI mask file (1 ROI per mask) with suffix, e.g., MyROIforFC+tlrc.  2nd column is the name (one argument) of the contrast you are using, e.g., Fear-Happy, if you want to examine a psychophysiological interaction.  If not, just name it 'NoContrast' or something like that.  The third column is the name of the seed ROI for the mask, e.g., InsulaSeed.  You can have multiple lines in this file if you want to examine multiple seed regions.
set FCtype = $11
#Set this to '1' if you are looking at a PPI (e.g., a contrast like Fear-Happy).  You will also need a 1 column file with 1's, -1's and 0's corresponding to the timepoints you wish to contrast or ignore (these should be the timepoints prior to convolving with the HRF).  Call this [task]FC[contrast].1D and put it in your scripts directory.  Also, you will need a task regressor (convolved with the HRF) to signify this contrast for the deconvolve step.  Call this [contrast].1D and put it in your scripts directory.  Set this to '2' if you want to look at connectivity across the entire timeseries without regard to task conditions.  If you want to do the latter (#2), you will need to modify the end of the script to enter a 3dDeconvolve command with the appropriate regressors for other conditions you are not interested in (e.g., the task regressors) in addition to the seed timeseries.
set censor = $12
#The name of the censor file you are using to censor outlier timepoints.  You only need to enter the suffix, not the subject's number (e.g., if censor file is 1111201201hariri_dtsali2anatCensor.1D, just enter hariri_dtsali2anatCensor.1D)
#
#
set group = $13
#group is either MDD or ctl
setenv FSLDIR /mnt/nfs/software/fsl-4.1.6
source ${FSLDIR}/etc/fslconf/fsl.csh
setenv PATH ${FSLDIR}/bin:${PATH}
#set study directory
set studydir = "${studydirectory}"
cd ${studydir}/scripts/
foreach subj (${subject})
cd ${studydir}/${group}/${subj}BRIKS/
if (${scripttype} == 1) then
rm ${subj}${task}_WMts+orig*
rm ${subj}${task}_LateralVentriclesTS.1D
rm ${subj}${SSsuffix}_VBMmask+orig*
rm ${subj}${SSsuffix}_FSLGMmask+orig*
rm ${subj}${SSsuffix}_FSLWMmask+orig*
rm ${subj}${SSsuffix}_FSLCSFmask+orig*
rm ${subj}${SSsuffix}_FSLGMmask4.0+orig*
rm ${subj}${SSsuffix}_FSLWMmask4.0+orig*
rm ${subj}${SSsuffix}_FSLCSFmask4.0+orig*
rm ${subj}${SSsuffix}_FSLWMmask4.0_eroded+orig*
rm ${subj}${SSsuffix}_FSLGMmask4.0Frac*+orig*
rm ${subj}${SSsuffix}_FSLCSFmask4.0AutoTal+tlrc*
rm ${subj}2ndVentricleMask_FSLCSFmask4.0AutoTal+tlrc*
rm ${subj}1stVentricleMask_FSLCSFmask4.0AutoTal+tlrc*
rm ${subj}LateralVentriclesMask_FSLCSFmask4.0AutoTal+tlrc*
rm ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+orig*
rm ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_eroded*+orig*
rm RPYXYZ.1D
rm ${subj}${task}_LHSinputs+tlrc*
rm ${subj}${task}_fitts+tlrc*
rm ${subj}${task}_residuals+tlrc*
rm ${subj}${task}_LHSinputs+orig*
rm ${subj}${task}_fitts+orig*
rm ${subj}${task}_residuals+orig*
rm *dewarped*
rm *xform*.1D
endif
#Check to see if talairached functional run exists
if (! -e ${subj}${task}+tlrc.HEAD) then
if (${taltype} == 1) then
@auto_tlrc -apar ${subj}${talsuffix}+tlrc -input ${subj}${task}+orig -dxyz 4.0 -rmode NN -verb
else
 adwarp -apar ${subj}${talsuffix}+tlrc \
	-dpar ${subj}${task}+orig\
	-prefix ${subj}${task}\
	-dxyz 4 \
	-force \
	-resam NN
endif
endif
#Set up file for FSL
if (! -e ${subj}${SSsuffix}anat.RPI.nii.gz) then
3dresample -orient RPI -inset ${studydir}/${group}/${subj}BRIKS/${subj}${SSsuffix}+orig -prefix ${subj}${SSsuffix}anat.RPI.nii
endif
#Run FAST in FSL
if (! -e ${subj}${SSsuffix}_VBMmask+orig.HEAD) then
/mnt/nfs/software/fsl-4.1.6/bin/fast -R 0.3 -H 0.1 ${subj}${SSsuffix}anat.RPI
#Turn output back into BRIK files
3dresample -orient ASL -inset ${subj}${SSsuffix}anat.RPI_seg.nii.gz -prefix ${subj}${SSsuffix}_VBMmask+orig
endif
#Create separate masks for each tissue type
if (! -e ${subj}${SSsuffix}_FSLGMmask+orig.HEAD) then
3dcalc -datum short -prefix ${subj}${SSsuffix}_FSLGMmask -a ${subj}${SSsuffix}_VBMmask+orig -expr "equals(a,2)"
3dcalc -datum short -prefix ${subj}${SSsuffix}_FSLWMmask -a ${subj}${SSsuffix}_VBMmask+orig -expr "equals(a,3)"
3dcalc -datum short -prefix ${subj}${SSsuffix}_FSLCSFmask -a ${subj}${SSsuffix}_VBMmask+orig -expr "equals(a,1)"
endif
#Resample each mask into 4x4x4 voxels
foreach mask (GMmask WMmask CSFmask)
if (! -e ${subj}${SSsuffix}_FSL${mask}4.0+orig.HEAD) then
3dresample -master ${subj}${task}+orig'['0']' -rmode NN -inset ${subj}${SSsuffix}_FSL${mask}+orig -prefix ${subj}${SSsuffix}_FSL${mask}4.0
endif
end
#Erode the WM mask
if (! -e ${subj}${SSsuffix}_FSLWMmask4.0_eroded+orig.HEAD) then
3dcalc -datum short -a ${subj}${SSsuffix}_FSLWMmask4.0+orig -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr "a*(1-amongst(0,b,c,d,e,f,g))" -prefix ${subj}${SSsuffix}_FSLWMmask4.0_eroded
endif
#Create optional fractionized GM masks
foreach clip (40 45 50)
if (! -e ${subj}${SSsuffix}_FSLGMmask4.0Frac${clip}+orig.HEAD) then
3dfractionize -template ${subj}${task}+orig -input ${subj}${SSsuffix}_FSLGMmask+orig -prefix ${subj}${SSsuffix}_FSLGMmask4.0Frac${clip} -clip ${clip} -vote
endif
end
#Create ventricle masks
if (! -e ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+orig.HEAD) then
if (${taltype} == 1) then
@auto_tlrc -apar ${subj}${talsuffix}+tlrc -input ${subj}${SSsuffix}_FSLCSFmask+orig -suffix 4.0AutoTal -dxyz 4.0 -rmode NN -verb
else
 adwarp -apar ${subj}${talsuffix}+tlrc \
	-dpar ${subj}${SSsuffix}_FSLCSFmask+orig\
	-prefix ${subj}${SSsuffix}_FSLCSFmask4.0AutoTal\
	-dxyz 4 \
	-force \
	-resam NN
endif
3dcalc -RAI -datum float -a ${subj}${SSsuffix}_FSLCSFmask4.0AutoTal+tlrc -expr "within(x,-28,0)*within(z,2,30)*within(y,-29,45)*step(a)" -prefix ${subj}2ndVentricleMask_FSLCSFmask4.0AutoTal
3dcalc -RAI -datum float -a ${subj}${SSsuffix}_FSLCSFmask4.0AutoTal+tlrc -expr "within(x,0,28)*within(z,2,30)*within(y,-29,45)*step(a)" -prefix ${subj}1stVentricleMask_FSLCSFmask4.0AutoTal
3dcalc -datum short -a ${subj}1stVentricleMask_FSLCSFmask4.0AutoTal+tlrc -b ${subj}2ndVentricleMask_FSLCSFmask4.0AutoTal+tlrc -expr "a+b" -prefix ${subj}LateralVentriclesMask_FSLCSFmask4.0AutoTal
cat_matvec ${subj}${talsuffix}+tlrc::WARP_DATA >> ${subj}${talsuffix}_xform.1D
3dWarp -matvec_out2in ${subj}${talsuffix}_xform.1D -NN -prefix ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+tlrc -gridset ${subj}${task}+orig ${subj}LateralVentriclesMask_FSLCSFmask4.0AutoTal+tlrc
3drefit -view orig ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+tlrc
endif
#Erode the ventricle masks (3 different types: 3 directions, 2 directions, 1 direction)
if (! -e ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx1d+orig.HEAD) then
3dcalc -datum short -a ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+orig -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k -expr "a*(1-amongst(0,b,c,d,e,f,g))" -prefix ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx3d
3dcalc -datum short -a ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+orig -b a+i -c a-i -d a+j -e a-j -expr "a*(1-amongst(0,b,c,d,e))" -prefix ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx2d
3dcalc -datum short -a ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0+orig -b a+i -c a-i -expr "a*(1-amongst(0,b,c))" -prefix ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx1d
endif
#Calculate the local WM regressors
if (! -e ${subj}${task}_WMts+orig.HEAD) then
3dLocalstat -nbhd 'SPHERE(30)' -stat mean -mask ${subj}${SSsuffix}_FSLWMmask4.0_eroded+orig -use_nonmask \
-prefix ${subj}${task}_WMts -datum float -quiet ${subj}${task}+orig
endif
#Calculate the ventricle timeseries
if (`3dinfo -VERB ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx1d+orig|grep "#"|awk '{print $12}'` != "0") then
set dim = 1
endif
if (`3dinfo -VERB ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx2d+orig|grep "#"|awk '{print $12}'` != "0") then
set dim = 2
endif
if (`3dinfo -VERB ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx3d+orig|grep "#"|awk '{print $12}'` != "0") then
set dim = 3
endif
if (! -e ${subj}${task}_LateralVentriclesTS.1D) then
3dROIstats -mask ${subj}${SSsuffix}_FSLLateralVentriclesmask4.0_erodedx${dim}d+orig -1Dformat -quiet ${subj}${task}+orig >> ${subj}${task}_LateralVentriclesTS.1D
endif
#Create the motion regre${SSsuffix}or timeseries
if (! -e RPYXYZ.1D) then
awk '{print $2,$3,$4,$5,$6,$7}' ${subj}${task}${regfile} >> RPYXYZ.1D
endif
cd ${studydir}/${group}/${subj}BRIKS/
if (${scripttype} == 1) then
rm ${subj}${task}_LHSinputs+tlrc*
rm ${subj}${task}_fitts+tlrc*
rm ${subj}${task}_residuals+tlrc*
rm ${subj}${task}_LHSinputs+orig*
rm ${subj}${task}_fitts+orig*
rm ${subj}${task}_residuals+orig*
rm RPYXYZ_det.1D
rm ${subj}${task}_WMts_det+orig*
rm ${subj}${task}_4thVentricleTS_det.1D
rm ${subj}${task}_3rdVentricleTS_det.1D
rm ${subj}${task}_2ndVentricleTS_det.1D
rm ${subj}${task}_1stVentricleTS_det.1D
rm ${subj}${task}_LateralVentriclesTS_det.1D
endif
#Detrend the motion, WM, and ventricle regressors
if (! -e RPYXYZ_det.1D) then
3dDetrend -DAFNI_1D_TRANOUT=YES -normalize -prefix RPYXYZ_det.1D -polort 1 RPYXYZ.1D\'
endif
if (! -e ${subj}${task}_WMts_det+orig.HEAD) then
3dDetrend -normalize -prefix ${subj}${task}_WMts_det -polort 1 ${subj}${task}_WMts+orig
endif
if (! -e ${subj}${task}_LateralVentriclesTS_det.1D) then
3dDetrend -DAFNI_1D_TRANOUT=YES -normalize -prefix ${subj}${task}_LateralVentriclesTS_det.1D -polort 1 ${subj}${task}_LateralVentriclesTS.1D\'
endif
#Create the fitted time series
if (! -e ${subj}${task}_LHSinputs+orig.HEAD) then
3dTfitter -RHS ${subj}${task}+orig -LHS RPYXYZ_det.1D ${subj}${task}_WMts_det+orig ${subj}${task}_LateralVentriclesTS_det.1D \
-prefix ${subj}${task}_LHSinputs -label roll pitch yaw x y z LocalWM LatVent -fitts ${subj}${task}_fitts -quiet
endif
if (${scripttype} == 1) then
rm ${subj}${task}_CleanedforFC+orig*
rm ${subj}${task}_CleanedforFC_*+tlrc*
rm ${subj}${task}_CleanedforFC_*+orig*
rm ${subj}${task}_CleanedforFC_FC*.1D
rm ${subj}${task}FC*
endif
#Create the residual time series
if (! -e ${subj}${task}_CleanedforFC+orig.HEAD) then
3dcalc -float -a ${subj}${task}+orig -b ${subj}${task}_fitts+orig -expr "a-b" -prefix ${subj}${task}_CleanedforFC+orig
endif
#Blur the residual time series
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts+tlrc.HEAD) then
3dBlurInMask \
-input ${subj}${task}_CleanedforFC+orig \
-FWHM ${blur} \
-mask ${subj}${SSsuffix}_FSLGMmask4.0Frac45+orig \
-preserve \
-prefix ${subj}${task}_CleanedforFC_GM${blur}blur \
-quiet
3dBlurInMask \
-input ${subj}${task}_CleanedforFC_GM${blur}blur+orig \
-FWHM ${blur} \
-mask ${subj}${SSsuffix}_FSLWMmask4.0_eroded+orig \
-preserve \
-prefix ${subj}${task}_CleanedforFC_GMWM${blur}blur \
-quiet
#Bandwidth filter the residual blurred timeseries
3dFourier \
-prefix ${subj}${task}_CleanedforFC_FC${blur}BFts \
-highpass .009 -lowpass .08 \
-retrend \
${subj}${task}_CleanedforFC_GMWM${blur}blur+orig
#Put the blurred, filtered, residual timeseries into talairach space
echo "The taltype is::::::::::::::::::::::::::::::::"
echo ${taltype}
if (${taltype} == 1) then
@auto_tlrc -apar ${subj}${talsuffix}+tlrc -input ${subj}${task}_CleanedforFC_FC${blur}BFts+orig -dxyz 4.0 -rmode NN -verb
else
 adwarp -apar ${subj}${talsuffix}+tlrc \
	-dpar ${subj}${task}_CleanedforFC_FC${blur}BFts+orig\
	-prefix ${subj}${task}_CleanedforFC_FC${blur}BFts\
	-dxyz 4 \
	-force \
	-resam NN
endif
endif
set reps = `3dnvals ${subj}${task}_CleanedforFC_FC${blur}BFts+tlrc`
foreach mask (`more ${studydir}/scripts/${roimasklist}|awk '{print $1}'`)
set contrast = `more ${studydir}/scripts/${roimasklist}|grep "${mask}"|awk '{print $2}'`
echo "Contrast:" 
echo ${contrast}
set tag = `more ${studydir}/scripts/${roimasklist}|grep "${mask}"|awk '{print $3}'`
cd ${studydir}/${group}/${subj}BRIKS/
echo "Tag:"
echo ${tag}
if (${scripttype} == 1) then
rm ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}*.1D
endif
echo "Herro"
#Extract the seed time series and create the PPI
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}.1D) then
3dROIstats -mask ${studydir}/testdir/${mask} -1Dformat -quiet ${subj}${task}_CleanedforFC_FC${blur}BFts+tlrc >> ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}.1D
endif
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detrend.1D) then
3dDetrend -polort 1 -prefix ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detrend.1D \
    ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}.1D\'
endif
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcol.1D) then
1dtranspose ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detrend.1D \
    ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcol.1D
endif
if (${FCtype} == 1) then
if (! -e ${studydir}/scripts/FChrf.1D) then
waver -TR 2.0 -GAM -inline 1@1 >> ${studydir}/scripts/FChrf.1D
endif
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolcon.1D) then
3dTfitter -RHS ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcol.1D \
    -FALTUNG ${studydir}/scripts/FChrf.1D ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolcon 012 -1
endif
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol.1D) then
1dtranspose ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolcon.1D \
    ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol.1D
endif
#Create the interaction term for the PPI
echo "HI THERE"
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_inter.1D) then
1deval -a ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol.1D \
    -b ${studydir}/scripts/${task}FC${contrast}.1D -expr "a*b" >> ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_inter.1D
endif
if (! -e ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_interwav.1D) then
waver -GAM -peak 1 -TR 2.0 -input ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_inter.1D -numout ${reps} >> ${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_interwav.1D
endif
endif
endif
if (${FCtype} == 1) then
if (! -e ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle.buc+orig.HEAD) then
3dDeconvolve \
-input ${subj}${task}_CleanedforFC_FC${blur}BFts+orig \
-censor ${subj}${censor} \
-num_stimts 3 \
-stim_file 1 ${studydir}/scripts/${contrast}.1D -stim_label 1 ${contrast} \
-stim_file 2 ${studydir}/${group}/${subj}BRIKS/${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcol.1D -stim_label 2 ROIts \
-stim_file 3 ${studydir}/${group}/${subj}BRIKS/${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcolconcol_interwav.1D \
-stim_label 3 ROIx${contrast} \
-rout -tout -bout -nofull_first -float \
-bucket ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle.buc
endif
#Fisher-Z transform the correlations; if manual taliaraching was done, will need to alter the talairaching command at the end of 3ddExtall2 to "adwarp" rather than "@auto_tlrc"
if (! -e ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyleallrFz+tlrc.HEAD) then
cd ${studydir}/scripts/
if (${taltype} == 1) then
echo "Doing fisher-Z correlations::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
./3ddExtallrFz ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir}/${group} 1 ${talsuffix}
./3ddExtallCoef ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir}/${group} 1 ${talsuffix}
else
./3ddExtallrFz ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir}/${group} 2 ${talsuffix}
./3ddExtallCoef ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir}/${group} 2 ${talsuffix}
endif
endif
endif
#Here you will need to modify the 3dDeconvolve command below to control for other regressors of no interest, e.g., task conditions, in addition to the seed timeseries.
if (${FCtype} == 2) then
if (! -e ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle.buc+orig.HEAD) then
3dDeconvolve \
-input ${subj}${task}_CleanedforFC_FC${blur}BFts+orig \
-censor ${subj}${censor} \
-num_stimts ? \
-stim_file 1 ${studydir}/${group}/${subj}BRIKS/${subj}${task}_CleanedforFC_FC${blur}BFts_${contrast}_${tag}detcol.1D -stim_label 1 ROIts \
-stim_file 2 ...
-rout -tout -bout -nofull_first -float \
-bucket ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle.buc
endif
#Fisher-Z transform the correlations; if manual taliaraching was done, will need to alter the talairaching command at the end of 3ddExtall2 to "adwarp" rather than "@auto_tlrc"
if (! -e ${subj}${task}FC${blur}BFts_${contrast}_${tag}_SPMstyleallrFz+tlrc.HEAD) then
cd ${studydir}/scripts/
if (${taltype} == 1) then
./3ddExtallrFz ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir} 1 ${talsuffix}
./3ddExtallCoef ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir} 1 ${talsuffix}
else
./3ddExtallrFz ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir} 2 ${talsuffix}
./3ddExtallCoef ${subj} ${task}FC${blur}BFts_${contrast}_${tag}_SPMstyle ${studydir} 2 ${talsuffix}
endif
endif
endif
end
end
