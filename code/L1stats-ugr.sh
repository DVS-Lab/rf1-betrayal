#!/usr/bin/env bash

# This script will perform Level 1 statistics in FSL.
# Rather than having multiple scripts, we are merging three analyses
# into this one script:
#		1) activation
#		2) seed-based ppi
#		3) network-based ppi
# Note that activation analysis must be performed first.
# Seed-based PPI should follow activation analyses.

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
rf1datadir=/ZPOOL/data/projects/rf1-sra-linux2 #need to fix this upon release (no hard coding paths)

# study-specific inputs
TASK=ugr
sm=5
sub=$1
run=$2
ppi=$3 # 0 for activation, otherwise seed region (check path on line 84; expectation: ${maindir}/masks/seed-${ppi}.nii.gz)
model=3 # 2 for the "original" merged events version
ses=`zeropad 1 2`

# set inputs and general outputs (may change depending on Tedana or fMRIPrep confounds)
MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}/ses-${ses}
mkdir -p $MAINOUTPUT
DATA=${rf1datadir}/derivatives/fmriprep-24/sub-${sub}/ses-${ses}/func/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
NVOLUMES=`fslnvols $DATA`
CONFOUNDEVS=${rf1datadir}/derivatives/fsl/confounds_tedana-24/sub-${sub}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv
# echo $CONFOUNDEVS
echo "Starting analysis for sub-${sub}, run-${run}, analysis type: ${ppi}" 
if [ ! -e $CONFOUNDEVS ]; then
	echo "missing confounds: $CONFOUNDEVS "  
	exit # exiting to ensure nothing gets run without confounds
fi

EVDIR=${maindir}/derivatives/fsl/EVfiles/sub-${sub}/${TASK}/model-3/run-${run}

# empty EVs (specific to this study)
# 10 indicates empty
MISSED_TRIAL=${EVDIR}_missed_trial.txt
if [ -e $MISSED_TRIAL ]; then
	SHAPE_EV=3
else
	SHAPE_EV=10
fi

# set output based in whether it is activation or ppi
if [ "$ppi" == "0" ]; then
        TYPE=act
        OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}_sm-${sm}
        OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${model}_type-${TYPE}_run-${run}.fsf
else
        TYPE=ppi
        OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_seed-${ppi}_run-${run}_sm-${sm}
        OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${model}_type-${TYPE}_seed-${ppi}_run-${run}.fsf
fi

# check for output and skip existing
if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
        exit
else
        rm -rf ${OUTPUT}.feat
fi

# create template and run analyses
ITEMPLATE=${maindir}/templates/L1_task-ugr_model-${model}_type-${TYPE}.fsf # this should be pointed to the trust template that Shenghan is using, which should be identical to the one I made long ago
if [ "$ppi" == "0" ]; then
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@DATA@'$DATA'@g' \
        -e 's@EVDIR@'$EVDIR'@g' \
        -e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
        -e 's@SHAPE_EV@'$SHAPE_EV'@g' \
        -e 's@SMOOTH@'$sm'@g' \
        -e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
        -e 's@NVOLUMES@'$NVOLUMES'@g' \
        <$ITEMPLATE> $OTEMPLATE
        feat $OTEMPLATE
else
        PHYS=${maindir}/derivatives/fsl/sub-${sub}/ses-${ses}/sub-${sub}_ses-${ses}_task-${TASK}_run-${run}_${ppi}.txt
        MASK=${maindir}/masks/seed-${ppi}.nii.gz
        fslmeants -i $DATA -o $PHYS -m $MASK
        sed -e 's@OUTPUT@'$OUTPUT'@g' \
        -e 's@DATA@'$DATA'@g' \
        -e 's@EVDIR@'$EVDIR'@g' \
        -e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
        -e 's@SHAPE_EV@'$SHAPE_EV'@g' \
        -e 's@PHYS@'$PHYS'@g' \
        -e 's@SMOOTH@'$sm'@g' \
        -e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
        -e 's@NVOLUMES@'$NVOLUMES'@g' \
        <$ITEMPLATE> $OTEMPLATE
        feat $OTEMPLATE
fi

## fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-preprocessed-by-fmriprep-without-re-registering-the-data-to-the-mni-space/784/3
mkdir -p ${OUTPUT}.feat/reg
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
cp $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
cp ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz


# delete unused files
rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz

	



