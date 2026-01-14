#!/bin/bash

BASE_DIR="/ZPOOL/data/projects/rf1-betrayal/derivatives/fsl"
OUTFILE="/ZPOOL/data/projects/rf1-betrayal/code/filepaths-trust.txt"

# clear output file if it exists
> "${OUTFILE}"

for sub in $(cat sublist_DD128.txt); do

  #paths for checking if the files exist
  L2_PATH_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L2_task-trust_model-01_type-act_sm-5.gfeat"
  L2_EVPATH_CHECK="${BASE_DIR}/EVfiles/sub-${sub}/trust/run-2_choice_computer.txt"
  L1_RUN1_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-1_sm-5.feat"
  L1_RUN2_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-2_sm-5.feat"


  L2_PATH="${BASE_DIR}/sub-${sub}/ses-01/L2_task-trust_model-01_type-act_sm-5.gfeat/copeCOPENUM.feat/stats/cope1.nii.gz"
  L1_RUN1="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-1_sm-5.feat/stats/copeCOPENUM.nii.gz"
  L1_RUN2="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-2_sm-5.feat/stats/copeCOPENUM.nii.gz"

	#echo ${L2_PATH}
	#echo ${L1_RUN2}
	#echo ${L1_PATH}
	
  if [ -d "${L1_RUN2_CHECK}" ] && [ -f "${L2_EVPATH_CHECK}" ]; then
    # run-2 exists → use L2
    [ -d "${L2_PATH_CHECK}" ] && echo "${L2_PATH}" >> "${OUTFILE}"

  else
    # run-2 does not exist → use L1 run-1
    [ -d "${L1_RUN1_CHECK}" ] && echo "${L1_RUN1}" >> "${OUTFILE}"
  fi

done


# if 'cat -A sublist_n241.txt | head' shows 10317^M$ then run 'sed -i 's/\r$//' sublist_n241.txt'
