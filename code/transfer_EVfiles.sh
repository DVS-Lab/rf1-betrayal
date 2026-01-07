#!/usr/bin/env bash

# source and destination base directories
SRC_BASE="/ZPOOL/data/projects/rf1-betrayal/dummy_EVfiles"
DST_BASE="/ZPOOL/data/projects/rf1-betrayal/derivatives/fsl/EVfiles"

# loop over subject directories in dummy_EVfiles
for subj_dir in "${SRC_BASE}"/sub-*; do
    subj=$(basename "${subj_dir}")

    SRC_UGR="${SRC_BASE}/${subj}/ugr"
    DST_SUBJ="${DST_BASE}/${subj}"

    # check that source ugr exists
    if [ ! -d "${SRC_UGR}" ]; then
        echo "Skipping ${subj}: no ugr directory in dummy_EVfiles"
        continue
    fi

    # check that destination subject directory exists
    if [ ! -d "${DST_SUBJ}" ]; then
        echo "Skipping ${subj}: destination subject directory does not exist"
        continue
    fi

    echo "Copying ugr for ${subj}"
    cp -r "${SRC_UGR}" "${DST_SUBJ}/"
done

echo "Done."
