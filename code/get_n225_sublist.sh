#!/bin/bash

# ===============================
# CONFIGURATION
# ===============================

rf1datadir=/ZPOOL/data/projects/rf1-sra-linux2

DERIV_ROOT=${rf1datadir}/derivatives/fmriprep-24
BEHAV_ROOT=/ZPOOL/data/projects/rf1-sra/stimuli

TASK1="trust"
TASK2="ugr"

SPACE="space-MNI152NLin6Asym"
DESC="desc-preproc_bold.nii.gz"

OUT_SUBJECT_LIST_ALL="subject_list_all.txt"
OUT_SUBJECT_LIST_FIRST225="subject_list_n225.txt"
OUT_EXCLUSION="exclusion_log.tsv"
OUT_QC="qc_table.tsv"

# ===============================
# INITIALIZE OUTPUT FILES
# ===============================

> "$OUT_SUBJECT_LIST_ALL"
echo -e "sub\treason" > "$OUT_EXCLUSION"
echo -e "sub\t${TASK1}_bold_ok\t${TASK2}_bold_ok\t${TASK1}_beh_ok\t${TASK2}_beh_ok\tinclude" > "$OUT_QC"

# ===============================
# LOOP THROUGH SUBJECTS
# ===============================

for sub_dir in ${DERIV_ROOT}/sub-*; do

    sub=$(basename "$sub_dir")
    sub_num=${sub#sub-}

    # Remove leading zeros (sub-001 → 1)
    sub_num=$(echo "$sub_num" | sed 's/^0*//')

    echo "Checking sub-${sub_num}"

    include=1
    reason=""

    trust_bold_ok=1
    ugr_bold_ok=1
    trust_beh_ok=1
    ugr_beh_ok=1

    # ---------------------------------------
    # TRUST neuroimaging check
    # ---------------------------------------

    trust_bolds=($(find "$sub_dir" -type f \
        -name "*task-${TASK1}_*${SPACE}*${DESC}" \
        ! -name "*echo-*" 2>/dev/null))

    if [ ${#trust_bolds[@]} -eq 0 ]; then
        trust_bold_ok=0
        include=0
        reason="${reason} Missing trust preproc bold;"
    fi

    # ---------------------------------------
    # UGR neuroimaging check
    # ---------------------------------------

    ugr_bolds=($(find "$sub_dir" -type f \
        -name "*task-${TASK2}_*${SPACE}*${DESC}" \
        ! -name "*echo-*" 2>/dev/null))

    if [ ${#ugr_bolds[@]} -eq 0 ]; then
        ugr_bold_ok=0
        include=0
        reason="${reason} Missing ugr preproc bold;"
    fi

    # ---------------------------------------
    # Behavioral checks
    # ---------------------------------------

    trust_dir="${BEHAV_ROOT}/Scan-Investment_Game/logs/${sub_num}"
    ugr_dir="${BEHAV_ROOT}/Scan-Lets_Make_A_Deal/logs/${sub_num}"

    trust_files=($(find "$trust_dir" -type f \
        -name "sub-${sub_num}_task-trust*_raw.csv" 2>/dev/null))

    if [ ${#trust_files[@]} -eq 0 ]; then
        trust_beh_ok=0
        include=0
        reason="${reason} Missing trust behavioral;"
    fi

    ugr_files=($(find "$ugr_dir" -type f \
        -name "sub-${sub_num}_task-ultimatum*_raw.csv" 2>/dev/null))

    if [ ${#ugr_files[@]} -eq 0 ]; then
        ugr_beh_ok=0
        include=0
        reason="${reason} Missing ugr behavioral;"
    fi

    # ---------------------------------------
    # Write QC table (numeric IDs)
    # ---------------------------------------

    echo -e "${sub_num}\t${trust_bold_ok}\t${ugr_bold_ok}\t${trust_beh_ok}\t${ugr_beh_ok}\t${include}" >> "$OUT_QC"

    # ---------------------------------------
    # Inclusion / exclusion
    # ---------------------------------------

    if [ "$include" -eq 1 ]; then
        echo "${sub_num}" >> "$OUT_SUBJECT_LIST_ALL"
    else
        echo -e "${sub_num}\t${reason}" >> "$OUT_EXCLUSION"
    fi

done

# ---------------------------------------
# Create first 225 subject list
# ---------------------------------------

sort -n "$OUT_SUBJECT_LIST_ALL" | head -n 225 > "$OUT_SUBJECT_LIST_FIRST225"

echo "Dual-task subject list generation complete."
echo "All valid subjects: $OUT_SUBJECT_LIST_ALL"
echo "First 225 subjects: $OUT_SUBJECT_LIST_FIRST225"
