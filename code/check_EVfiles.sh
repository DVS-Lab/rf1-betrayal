#!/bin/bash

BASE_DIR="/ZPOOL/data/projects/rf1-betrayal/derivatives/fsl"
SUBLIST="subject_list_n225.txt"
OUTCSV="/ZPOOL/data/projects/rf1-betrayal/code/ev_qc_dualtask.csv"

# Write header
echo "sub_id,trust_run1_ok,trust_run2_ok,trust_valid_runs,ugr_run1_ok,ugr_run2_ok,ugr_valid_runs" > "${OUTCSV}"

for sub in $(cat "$SUBLIST"); do

  # ===============================
  # TRUST EV CHECKS
  # ===============================

  TRUST_RUN1="${BASE_DIR}/EVfiles/sub-${sub}/trust/run-1_choice_computer.txt"
  TRUST_RUN2="${BASE_DIR}/EVfiles/sub-${sub}/trust/run-2_choice_computer.txt"

  trust_run1_ok=0
  trust_run2_ok=0

  [ -s "$TRUST_RUN1" ] && trust_run1_ok=1
  [ -s "$TRUST_RUN2" ] && trust_run2_ok=1

  if [ "$trust_run1_ok" -eq 1 ] && [ "$trust_run2_ok" -eq 1 ]; then
      trust_status="both"
  elif [ "$trust_run1_ok" -eq 1 ]; then
      trust_status="run1"
  elif [ "$trust_run2_ok" -eq 1 ]; then
      trust_status="run2"
  else
      trust_status="none"
  fi

  # ===============================
  # UGR EV CHECKS
  # ===============================

  UGR_RUN1="${BASE_DIR}/EVfiles/sub-${sub}/ugr/model-3/run-1_nonsocial_high_pmod.txt"
  UGR_RUN2="${BASE_DIR}/EVfiles/sub-${sub}/ugr/model-3/run-2_nonsocial_high_pmod.txt"

  ugr_run1_ok=0
  ugr_run2_ok=0

  [ -s "$UGR_RUN1" ] && ugr_run1_ok=1
  [ -s "$UGR_RUN2" ] && ugr_run2_ok=1

  if [ "$ugr_run1_ok" -eq 1 ] && [ "$ugr_run2_ok" -eq 1 ]; then
      ugr_status="both"
  elif [ "$ugr_run1_ok" -eq 1 ]; then
      ugr_status="run1"
  elif [ "$ugr_run2_ok" -eq 1 ]; then
      ugr_status="run2"
  else
      ugr_status="none"
  fi

  # ===============================
  # WRITE OUTPUT
  # ===============================

  echo "${sub},${trust_run1_ok},${trust_run2_ok},${trust_status},${ugr_run1_ok},${ugr_run2_ok},${ugr_status}" >> "${OUTCSV}"

done

echo "Dual-task EV QC written to ${OUTCSV}"
