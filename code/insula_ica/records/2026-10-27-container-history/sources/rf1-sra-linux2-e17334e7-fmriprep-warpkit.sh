#!/usr/bin/env bash

# Sequential OpenNeuro smoke tests for the fMRIPrep warpkit/MEDIC branch.
# The script installs each dataset with DataLad, downloads one participant,
# runs the local Apptainer image, and then moves to the next dataset.

set -euo pipefail

# Ensure paths are correct irrespective of where the user runs the script.
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
maindir="$(dirname "$scriptdir")"

# Dedicated validation sandbox. All downloaded inputs, outputs, logs, and work
# files stay under this folder so the existing project tree remains uncluttered.
testrunroot="$maindir/derivatives/fmriprep-warpkit-validation"
bidsroot="$testrunroot/bids"
outputroot="$testrunroot/outputs"
logroot="$testrunroot/logs"
scratchroot="$testrunroot/work"
templateflow_host=/ZPOOL/data/tools/templateflow
mplconfig_host=/ZPOOL/data/tools/mplconfigdir
licenses_host=/ZPOOL/data/tools/licenses
fmriprep_img=/ZPOOL/data/tools/fmriprep-warpkit-andrew.sif # made from https://github.com/DVSneuro/fmriprep/tree/warpkit-medic

# Container paths.
base_container=/base
templateflow_container=/opt/templateflow
mplconfig_container=/opt/mplconfigdir
licenses_container=/opts

# OpenNeuro datasets with complex-valued multi-echo data suggested for testing.
# Your datasets are intentionally omitted here: ds005123, ds006707, ds007486.
# Format: dataset_id version_tag subject_selector
# The selector may be a full subject label or a 1-based index into sorted sub-* directories.
# ds006072/sub-P1 fails appropriately because it has no matching part-phase BOLD files.

#datasets=(
#    "ds007637 1.0.0 sub-01"
#    "ds006926 1.0.0 sub-a01"
#    "ds006131 2.0.1 sub-20188"
#    "ds006072 1.0.8 sub-P1R"
#    "ds005085 1.0.0 sub-10015"
#    "ds002278 2.0.0 sub-Blossom"
#)

datasets=(
    "ds006072 1.0.8 sub-P1R"
)

run_id="$(date +%Y%m%d-%H%M%S)"

mkdir -p "$bidsroot" "$outputroot" "$logroot" "$scratchroot"

export APPTAINERENV_TEMPLATEFLOW_HOME=$templateflow_container
export APPTAINERENV_MPLCONFIGDIR=$mplconfig_container

install_dataset() {
    local dataset_id=$1
    local dataset_version=$2
    local dataset_dir="$bidsroot/$dataset_id"
    local source_url="https://github.com/OpenNeuroDatasets/${dataset_id}.git"

    if [[ -d "$dataset_dir/.datalad" ]]; then
        echo "[$dataset_id] Dataset already installed at $dataset_dir"
        if [[ "$dataset_version" == "latest" ]]; then
            datalad -C "$dataset_dir" update --merge
        else
            git -C "$dataset_dir" fetch --tags origin
        fi
    else
        echo "[$dataset_id] Installing from $source_url"
        datalad install -s "$source_url" "$dataset_dir"
        if [[ "$dataset_version" != "latest" ]]; then
            git -C "$dataset_dir" fetch --tags origin
        fi
    fi

    if [[ "$dataset_version" != "latest" ]]; then
        echo "[$dataset_id] Checking out OpenNeuro version $dataset_version"
        git -C "$dataset_dir" checkout --quiet "$dataset_version"
    fi
}

select_subject() {
    local dataset_id=$1
    local subject_selector=$2
    local dataset_dir="$bidsroot/$dataset_id"
    local subject_index
    local subjects=()

    if [[ "$subject_selector" == sub-* ]]; then
        if [[ -d "$dataset_dir/$subject_selector" ]]; then
            printf '%s\n' "$subject_selector"
            return 0
        fi

        echo "[$dataset_id] Could not find requested subject $subject_selector" >&2
        return 1
    fi

    subject_index=$((subject_selector - 1))
    mapfile -t subjects < <(find "$dataset_dir" -maxdepth 1 -type d -name 'sub-*' -exec basename {} \; | sort)
    if (( ${#subjects[@]} <= subject_index )); then
        echo "[$dataset_id] Could not find sorted subject number $subject_selector" >&2
        return 1
    fi

    printf '%s\n' "${subjects[$subject_index]}"
}

get_subject_data() {
    local dataset_id=$1
    local subject=$2
    local dataset_dir="$bidsroot/$dataset_id"

    echo "[$dataset_id] Downloading ${subject}"
    datalad -C "$dataset_dir" get dataset_description.json || true
    datalad -C "$dataset_dir" get participants.tsv participants.json || true
    datalad -C "$dataset_dir" get -r "$subject"
}

run_fmriprep() {
    local dataset_id=$1
    local dataset_version=$2
    local subject=$3
    local participant_label="${subject#sub-}"
    local out_dir="$outputroot/$dataset_id"
    local work_dir="$scratchroot/$dataset_id"
    local version_label="${dataset_version//./_}"
    local log_file="$logroot/${dataset_id}_v${version_label}_${subject}_run-${run_id}_fmriprep-warpkit.log"

    mkdir -p "$out_dir" "$work_dir"

    echo "[$dataset_id] Running fMRIPrep for ${subject} from OpenNeuro version ${dataset_version}"
    echo "[$dataset_id] Log: $log_file"

    apptainer run --cleanenv \
        -B "${templateflow_host}:${templateflow_container}" \
        -B "${mplconfig_host}:${mplconfig_container}" \
        -B "${testrunroot}:${base_container}" \
        -B "${scriptdir}:${base_container}/code" \
        -B "${licenses_host}:${licenses_container}" \
        "${fmriprep_img}" \
        "${base_container}/bids/${dataset_id}" "${base_container}/outputs/${dataset_id}" \
        participant --participant_label "$participant_label" \
        --session-label 1 \
        --stop-on-first-crash \
        --me-output-echos \
        --output-spaces MNI152NLin6Asym \
        --bids-filter-file "${base_container}/code/fmriprep_config.json" \
        --skip-bids-validation \
        --fs-no-reconall \
        --fs-license-file "${licenses_container}/fs_license.txt" \
        --me-use-warpkit \
        --me-warpkit-noise-frames 3 \
        --nprocs 16 \
        --omp-nthreads 6 \
        --mem-mb 140000 \
        --low-mem \
        --resource-monitor \
        -w "${base_container}/work/${dataset_id}" \
        2>&1 | tee "$log_file"
}

failures=()

for dataset_spec in "${datasets[@]}"; do
    read -r dataset_id dataset_version subject_selector <<< "$dataset_spec"

    echo
    echo "========== $dataset_id $dataset_version =========="

    if ! install_dataset "$dataset_id" "$dataset_version"; then
        echo "[$dataset_id] DataLad install/update failed" >&2
        failures+=("$dataset_id: install")
        continue
    fi

    if ! subject="$(select_subject "$dataset_id" "$subject_selector")"; then
        failures+=("$dataset_id: subject selection")
        continue
    fi

    if ! get_subject_data "$dataset_id" "$subject"; then
        echo "[$dataset_id] DataLad get failed for $subject" >&2
        failures+=("$dataset_id/$subject: datalad get")
        continue
    fi

    if ! run_fmriprep "$dataset_id" "$dataset_version" "$subject"; then
        echo "[$dataset_id] fMRIPrep failed for $subject" >&2
        failures+=("$dataset_id/$subject: fmriprep")
        continue
    fi

    echo "[$dataset_id] Completed $subject"
done

echo
echo "========== Summary =========="
if (( ${#failures[@]} == 0 )); then
    echo "All dataset runs completed successfully."
else
    echo "Some dataset runs failed:"
    printf '  - %s\n' "${failures[@]}"
    exit 1
fi
