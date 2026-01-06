#!/usr/bin/env bash

# This script will perform Level 1 statistics in FSL.
# Rather than having multiple scripts, we are merging three analyses
# into this one script:
#		1) activation
#		2) seed-based ppi
#		3) network-based ppi
# Note that activation analysis must be performed first.
# Seed-based PPI and Network PPI should follow activation analyses.

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"
rf1datadir=/ZPOOL/data/projects/rf1-sra-data #need to fix this upon release (no hard coding paths)

# study-specific inputs
TASK=ugr
sm=5
sub=$1
run=$2
ppi=$3 # 0 for activation, otherwise seed region or network
model=3 # 2 for the "original" merged events version
# maskname=$4

# set inputs and general outputs (may change depending on Tedana or fMRIPrep confounds)
MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}
mkdir -p $MAINOUTPUT
DATA=${rf1datadir}/derivatives/fmriprep/sub-${sub}/func/sub-${sub}_task-${TASK}_run-${run}_part-mag_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz
NVOLUMES=`fslnvols $DATA`
CONFOUNDEVS=${rf1datadir}/derivatives/fsl/confounds_tedana/sub-${sub}/sub-${sub}_task-${TASK}_run-${run}_desc-TedanaPlusConfounds.tsv
# echo $CONFOUNDEVS
echo "Starting analysis for sub-${sub}, run-${run}, analysis type: ${ppi}" 
if [ ! -e $CONFOUNDEVS ]; then
	echo "missing confounds: $CONFOUNDEVS "  
	exit # exiting to ensure nothing gets run without confounds
fi

EVDIR=${maindir}/derivatives/fsl/EVfiles/sub-${sub}/ugr/model-${model}/run-${run}

# empty EVs (specific to this study)
# 10 indicates empty
MISSED_TRIAL=${EVDIR}_missed_trial.txt
if [ -e $MISSED_TRIAL ]; then
	SHAPE_EV=3
else
	SHAPE_EV=10
fi


# if network (ecn or dmn), do nppi; otherwise, do activation or seed-based ppi
if [ "$ppi" == "ecn" -o  "$ppi" == "dmn" ]; then

	# check for output and skip existing
	OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-nppi-${ppi}_run-${run}_sm-${sm}
	if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
		echo "skipping sub-${sub}"
		exit
	else
	#	echo "running: $OUTPUT " 
		rm -rf ${OUTPUT}.feat
	fi

	# network extraction. need to ensure you have run Level 1 activation
	MASK=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-act_run-${run}_sm-${sm}.feat/mask
	if [ ! -e ${MASK}.nii.gz ]; then
		echo "cannot run nPPI because you're missing $MASK"
		exit
	fi
	for net in `seq 0 9`; do
		NET=${maindir}/masks/nan_rPNAS_2mm_net000${net}.nii.gz
		TSFILE=${MAINOUTPUT}/ts_task-${TASK}_net000${net}_nppi-${ppi}_run-${run}.txt
		fsl_glm -i $DATA -d $NET -o $TSFILE --demean -m $MASK
		eval INPUT${net}=$TSFILE
	done

	# set names for network ppi (we generally only care about ECN and DMN)
	DMN=$INPUT3
	ECN=$INPUT7
	if [ "$ppi" == "dmn" ]; then
		MAINNET=$DMN
		OTHERNET=$ECN
	else
		MAINNET=$ECN
		OTHERNET=$DMN

	# create template and run analyses
		ITEMPLATE=${maindir}/templates/L1_task-${TASK}_model-${model}_type-nppi.fsf
		OTEMPLATE=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_seed-${ppi}_run-${run}.fsf
			sed -e 's@OUTPUT@'$OUTPUT'@g' \
			-e 's@DATA@'$DATA'@g' \
			-e 's@EVDIR@'$EVDIR'@g' \
			-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
			-e 's@SHAPE_EV@'$SHAPE_EV'@g' \
			-e 's@EV_FRIENDN@'$EV_FRIENDN'@g' \
			-e 's@SHAPE_FRIENDN@'$SHAPE_FRIENDN'@g' \
			-e 's@EV_COMPN@'$EV_COMPN'@g' \
			-e 's@SHAPE_COMPN@'$SHAPE_COMPN'@g' \
			-e 's@EV_STRANGERN@'$EV_STRANGERN'@g' \
			-e 's@SHAPE_STRANGERN@'$SHAPE_STRANGERN'@g' \
			-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
			-e 's@MAINNET@'$MAINNET'@g' \
			-e 's@OTHERNET@'$OTHERNET'@g' \
			-e 's@INPUT0@'$INPUT0'@g' \
			-e 's@INPUT1@'$INPUT1'@g' \
			-e 's@INPUT2@'$INPUT2'@g' \
			-e 's@INPUT4@'$INPUT4'@g' \
			-e 's@INPUT5@'$INPUT5'@g' \
			-e 's@INPUT6@'$INPUT6'@g' \
			-e 's@INPUT8@'$INPUT8'@g' \
			-e 's@INPUT9@'$INPUT9'@g' \
			-e 's@INPUT10@'$INPUT10'@g' \
			-e 's@INPUT11@'$INPUT11'@g' \
			-e 's@NVOLUMES@'$NVOLUMES'@g' \
			<$ITEMPLATE> $OTEMPLATE
			feat $OTEMPLATE
	fi

else # otherwise, do activation and seed-based ppi

	# set output based in whether it is activation or ppi
	if [ "$ppi" == "0" ]; then
		TYPE=act
		OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_run-${run}_sm-${sm}
		OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${model}_type-${TYPE}_run-${run}.fsf
		name=act
	else
		TYPE=ppi
		OUTPUT=${MAINOUTPUT}/L1_task-${TASK}_model-${model}_type-${TYPE}_seed-${ppi}_run-${run}_sm-${sm}
		OTEMPLATE=${MAINOUTPUT}/L1_sub-${sub}_task-${TASK}_model-${model}_type-${TYPE}_seed-${ppi}_run-${run}.fsf
		echo "Output dir: $OUTPUT"
		type=seed-${ppi}
	fi

	# check for output and skip existing
	if [ -e ${OUTPUT}.feat/cluster_mask_zstat1.nii.gz ]; then
		exit
	else
		echo "running: $OUTPUT " 
		rm -rf ${OUTPUT}.feat
	fi

	# create template and run analyses
	ITEMPLATE=${maindir}/templates/L1_task-ugr_model-${model}_type-${TYPE}.fsf
	echo $ITEMPLATE
	
	echo $OTEMPLATE
	if [ "$ppi" == "0" ]; then
		sed -e 's@OUTPUT@'$OUTPUT'@g' \
		-e 's@DATA@'$DATA'@g' \
		-e 's@EVDIR@'$EVDIR'@g' \
		-e 's@MISSED_TRIAL@'$MISSED_TRIAL'@g' \
		-e 's@SHAPE_EV@'$SHAPE_EV'@g' \
		-e 's@EV_FRIENDN@'$EV_FRIENDN'@g' \
		-e 's@SHAPE_FRIENDN@'$SHAPE_FRIENDN'@g' \
		-e 's@EV_COMPN@'$EV_COMPN'@g' \
		-e 's@SHAPE_COMPN@'$SHAPE_COMPN'@g' \
		-e 's@EV_STRANGERN@'$EV_STRANGERN'@g' \
		-e 's@SHAPE_STRANGERN@'$SHAPE_STRANGERN'@g' \
		-e 's@SMOOTH@'$sm'@g' \
		-e 's@CONFOUNDEVS@'$CONFOUNDEVS'@g' \
		-e 's@NVOLUMES@'$NVOLUMES'@g' \
		<$ITEMPLATE> $OTEMPLATE
	else
		PHYS=${maindir}/derivatives/fsl/sub-${sub}/sub-${sub}_task-${TASK}_run-${run}_${ppi}.txt
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
fi
# fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-preprocessed-by-fmriprep-without-re-registering-the-data-to-the-mni-space/784/3
#mkdir -p ${OUTPUT}.feat/reg
#ln -s $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
#ln -s $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
#ln -s ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz



## fix registration as per NeuroStars post:
# https://neurostars.org/t/performing-full-glm-analysis-with-fsl-on-the-bold-images-preprocessed-by-fmriprep-without-re-registering-the-data-to-the-mni-space/784/3
mkdir -p ${OUTPUT}.feat/reg
ln -s $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/example_func2standard.mat
ln -s $FSLDIR/etc/flirtsch/ident.mat ${OUTPUT}.feat/reg/standard2example_func.mat
ln -s ${OUTPUT}.feat/mean_func.nii.gz ${OUTPUT}.feat/reg/standard.nii.gz
#fi


# delete unused files
rm -rf ${OUTPUT}.feat/stats/res4d.nii.gz
rm -rf ${OUTPUT}.feat/stats/corrections.nii.gz
rm -rf ${OUTPUT}.feat/stats/threshac1.nii.gz
rm -rf ${OUTPUT}.feat/filtered_func_data.nii.gz

