#!/bin/bash

# example code for FMRIPREP
# runs FMRIPREP on input subject
# usage: bash run_fmriprep.sh sub
# example: bash run_fmriprep.sh 102

sub=$1

# ensure paths are correct irrespective of where the user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# host-system paths
bidsdir="$maindir/bids"
derivdir="$maindir/derivatives"
scratchdir=/ZPOOL/data/scratch/$(whoami)
templateflow_host=/ZPOOL/data/tools/templateflow
mplconfig_host=/ZPOOL/data/tools/mplconfigdir
licenses_host=/ZPOOL/data/tools/licenses
fmriprep_img=/ZPOOL/data/tools/fmriprep-warpkit.sif # made from a special branch on my fork: https://github.com/DVSneuro/fmriprep/tree/warpkit-medic

# container paths
base_container=/base
templateflow_container=/opt/templateflow
mplconfig_container=/opt/mplconfigdir
licenses_container=/opts
scratch_container=/scratch

# make derivatives folder if it doesn't exist
if [ ! -d "$derivdir" ]; then
    mkdir -p "$derivdir"
fi

# make scratch folder if it doesn't exist
if [ ! -d "$scratchdir" ]; then
    mkdir -p "$scratchdir"
fi

export APPTAINERENV_TEMPLATEFLOW_HOME=$templateflow_container
export APPTAINERENV_MPLCONFIGDIR=$mplconfig_container

apptainer run --cleanenv \
    -B ${templateflow_host}:${templateflow_container} \
    -B ${mplconfig_host}:${mplconfig_container} \
    -B ${maindir}:${base_container} \
    -B ${licenses_host}:${licenses_container} \
    -B ${scratchdir}:${scratch_container} \
    ${fmriprep_img} \
    ${base_container}/bids ${base_container}/derivatives/fmriprep-warpkit-testing \
    participant --participant_label $sub \
    --stop-on-first-crash \
    --me-output-echos \
    --output-spaces MNI152NLin6Asym \
    --bids-filter-file ${base_container}/code/fmriprep_config.json \
    --skip-bids-validation \
    --fs-no-reconall \
    --fs-license-file ${licenses_container}/fs_license.txt \
    --me-use-warpkit \
    -w ${scratch_container}
