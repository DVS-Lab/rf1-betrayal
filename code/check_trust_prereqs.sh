#!/bin/bash

BASE_DIR="/ZPOOL/data/projects/rf1-betrayal/derivatives/fsl"
OUTCSV="/ZPOOL/data/projects/rf1-betrayal/code/trust_model_qc.csv"

# write header
echo "sub_id,trust_valid_runs" > "${OUTCSV}"

for sub in $(cat subject_list_n225.txt); do

  # --------- CHECK FILES ----------
  L2_PATH_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L2_task-trust_model-01_type-act_sm-5.gfeat/cope1.feat/stats/cope1.nii.gz"
  EVPATH_CHECK1="${BASE_DIR}/EVfiles/sub-${sub}/trust/run-1_choice_computer.txt"
  EVPATH_CHECK2="${BASE_DIR}/EVfiles/sub-${sub}/trust/run-2_choice_computer.txt"
  L1_RUN1_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-1_sm-5.feat/stats/cope1.nii.gz"
  L1_RUN2_CHECK="${BASE_DIR}/sub-${sub}/ses-01/L1_task-trust_model-01_type-act_run-2_sm-5.feat/stats/cope1.nii.gz"

echo "Checking: ${EVPATH_CHECK2}"

  # --------- DECISION LOGIC ----------
  if [ -f "${L2_PATH_CHECK}" ] && \
     [ -f "${EVPATH_CHECK1}" ] && \
     [ -f "${EVPATH_CHECK2}" ]; then

      MODEL="both"

  elif [ -f "${L1_RUN1_CHECK}" ] && \
       [ -f "${EVPATH_CHECK1}" ]; then

      MODEL="run1"

  elif [ -f "${L1_RUN2_CHECK}" ] && \
       [ -f "${EVPATH_CHECK2}" ]; then

      MODEL="run2"

  else
      MODEL="none"
  fi

  echo "${sub},${MODEL}" >> "${OUTCSV}"

done

echo "Wrote QC results to ${OUTCSV}"
