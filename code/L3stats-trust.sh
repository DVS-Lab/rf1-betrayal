#!/bin/bash

# This script will perform Level 3 statistics in FSL.
# Rather than having multiple scripts, we are merging three analyses
# into this one script:
#		1) two groups (older vs. younger)
#		2) two groups (older vs. younger), with covariates
#		3) single group average
#
# This script can also run randomise (permutation-based stats) on existing output.
# By default, randomise will not be be run if FEAT analyses do not exist. In addition,
# randomise will only be run on copes above a specified number (see copenum_thresh_randomise variable).
# If you have no intention of running randomise, you set copenum_thresh_randomise=20 (> max of 19 copes)
# and you could uncomment out the rm lines that remove the filtered_func_data file (save disk space).

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# study-specific inputs and general output folder

copenum=$1
copename=$2
analysis=$3
REPLACEME=act # act  ppi, nppi # this defines the parts of the path that differ across analyses
type=${REPLACEME} # For output template

# Variables that change per analysis. Check carefully! 
#covariate=$4
# covariate=dACC
covariate=ONES
# Covariates:
# dACC

# ONESUNDER55
# ONES

# Alpha-RelDep-Int-wCovars 

# Alpha-ExtBet-Int-wCovars 
# Alpha-ExtWor-Int-wCovars 
# Alpha-IntBet-Int-wCovars 
# Alpha-IntWor-Int-wCovars 

# RelDep-ExtBet-Int-wCovars D
# RelDep-ExtWor-Int-wCovars D
# RelDep-IntBet-Int-wCovars  D
# RelDep-IntWor-Int-wCovars D
# USI-ADI-ExtBet-wCovars D

# USI-ADI-ExtBet-Int-wCovars n


# RelDep-ExtBet-MainEffects-wCovars change to 112 
# RelDep-ExtWor-MainEffects-wCovars 
# RelDep-IntBet-MainEffects-wCovars  
# RelDep-IntWor-MainEffects-wCovars 


# eventually delete
# Alpha-ExtWor-RelDep-Int-wCovars
# Alpha-IntBet-RelDep-Int-wCovars
# Alpha-IntWor-RelDep-Int-wCovars

N=209 # update with total n after exclusions

if [[ $analysis == "act" ]]; then
	template=L3_task-ugr_model-3_type-act_group-${covariate}_n${N}_flame1.fsf #
else
	template=L3_task-ugr_model-3_type-ppi_group-${covariate}_n${N}_flame1.fsf
fi

# Templates:
# L3_task-ugr_group_dACC_n118_flame1.fsf 
# L3_task-ugr_type-act_group-RelDep_n94_flame1.fsf (no model in name = model-2, added model-3 for most recent analysis)
# L3_task-ugr_type-act_group-ExtBet-ADI-USI_n94_flame1.fsf
# L3_task-ugr_type-act_group-ExtBet-RelDep-Int_n94_flame1.fsf

# Set once and then forget.
model=3 # or model 2
task=ugr
modeltype=flame1
templatedir="/ZPOOL/data/projects/rf1-norms/templates"
MAINOUTPUT=${maindir}/derivatives/fsl/L3-act/L3_model-${model}_task-${task}_type-${type}-n${N}-cov-${covariate}-${modeltype}

mkdir -p $MAINOUTPUT

# set outputs and check for existing
cnum_pad=`zeropad ${copenum} 2`
OUTPUT=${MAINOUTPUT}/L3_model-${model}_task-${task}_n${N}-cov-${covariate}-cope-${copenum}_cname-${copename}-${modeltype}

	echo "[$(date)] re-doing: ${OUTPUT}" >> re-runL3.log
	rm -rf ${OUTPUT}.gfeat

	# create template and run FEAT analyses
	ITEMPLATE=${templatedir}/${template}
	OTEMPLATE=${MAINOUTPUT}/L3_task-${task}_model-${model}_type-${type}_cope-${copenum}_cname-${copename}_${covariate}_n${N}_${modeltype}.fsf
	sed -e 's@OUTPUT@'$OUTPUT'@g' \
	-e 's@COPENUM@'$copenum'@g' \
	-e 's@REPLACEME@'$REPLACEME'@g' \
	-e 's@BASEDIR@'$maindir'@g' \
	<$ITEMPLATE> $OTEMPLATE
	/usr/local/fsl/bin/feat $OTEMPLATE
	# feat $OTEMPLATE

# delete unused files
 rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/res4d.nii.gz
 rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/corrections.nii.gz
 rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/threshac1.nii.gz
#rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/filtered_func_data.nii.gz
# rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/var_filtered_func_data.nii.gz
