#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# setting inputs and common variables
sub=$1
type=$2 # change to "act" ## EDIT THESE FOR MASK "ppi- " "nppi-"
task=ugr 
sm=5 # edit if necessary, smoothing kernel, should not change 
MAINOUTPUT=${maindir}/derivatives/fsl/sub-${sub}
model=3 # model 2 is "old" version (merged events), 3 is new (timing fixes)

NCOPES=17 # our model has 17 copes
# ppi has more contrasts than act (phys), so need a different L2 template?
if [ "${type}" == "act" ]; then
	ITEMPLATE=${maindir}/templates/L2_task-${task}_model-${model}_type-act.fsf
	NCOPES=${NCOPES}
else
	ITEMPLATE=${maindir}/templates/L2_task-${task}_model-${model}_type-ppi.fsf
	let NCOPES=${NCOPES}+1 # add 1 since we tend to only have one extra contrast for PPI (for PHYS)
fi

INPUT1=${MAINOUTPUT}/L1_task-${task}_model-${model}_type-${type}_run-1_sm-${sm}.feat
INPUT2=${MAINOUTPUT}/L1_task-${task}_model-${model}_type-${type}_run-2_sm-${sm}.feat

# check for existing output and re-do if missing/incomplete
OUTPUT=${MAINOUTPUT}/L2_task-${task}_model-${model}_type-${type}_sm-${sm}

if [ -e ${OUTPUT}.gfeat/cope${NCOPES}.feat/cluster_mask_zstat1.nii.gz ]; then # check last (act) or penultimate (ppi) cope
	echo "skipping existing output"
else
	echo "re-doing: ${OUTPUT}" >> re-runL2.log
	rm -rf ${OUTPUT}.gfeat

	# set output template and run template-specific analyses
	OTEMPLATE=${MAINOUTPUT}/L2_task-${task}_model-${model}_type-${type}.fsf
	sed -e 's@OUTPUT@'$OUTPUT'@g' \
	-e 's@INPUT1@'$INPUT1'@g' \
	-e 's@INPUT2@'$INPUT2'@g' \
	<$ITEMPLATE> $OTEMPLATE
	feat $OTEMPLATE

	# delete unused files
	for cope in `seq ${NCOPES}`; do
		rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/res4d.nii.gz
		rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/corrections.nii.gz
		rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/stats/threshac1.nii.gz
		# rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/filtered_func_data.nii.gz
		# rm -rf ${OUTPUT}.gfeat/cope${cope}.feat/var_filtered_func_data.nii.gz
	done

fi
