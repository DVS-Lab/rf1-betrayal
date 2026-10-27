#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
Usage: bash fmriprep.sh [--dry-run] [--overwrite] SUBJECT
USAGE
}

scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
# shellcheck source=code/pipeline_common.sh
source "${scriptdir}/pipeline_common.sh"
rf1_load_config

dry_run=0
overwrite=0
while (($#)); do
  case "$1" in
    --dry-run)
      dry_run=1
      shift
      ;;
    --overwrite)
      overwrite=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      usage
      exit 2
      ;;
    *)
      break
      ;;
  esac
done

if (($# != 1)); then
  usage
  exit 2
fi

sub="$1"
bidsdir="${PROJECT_ROOT}/bids"
derivdir="${PROJECT_ROOT}/derivatives"
freesurferdir="${derivdir}/freesurfer"
scratchdir="${SCRATCH_ROOT}/$(whoami)"
IFS=' ' read -r -a output_spaces <<< "$FMRIPREP_OUTPUT_SPACES"

rf1_require_dir "${bidsdir}/sub-${sub}"

if [[ "$overwrite" -ne 1 ]] && python3 "${scriptdir}/check_pipeline_state.py" fmriprep-complete "$bidsdir" "$derivdir" "$sub" >/dev/null; then
  echo "sub-${sub} already has practical fMRIPrep completion outputs; skipping"
  exit 0
fi

export APPTAINERENV_TEMPLATEFLOW_HOME=/opt/templateflow
export APPTAINERENV_MPLCONFIGDIR=/opt/mplconfigdir

cmd=(
  singularity run --cleanenv
  -B "${TEMPLATEFLOW_HOME}:/opt/templateflow"
  -B "${MPLCONFIGDIR_HOST}:/opt/mplconfigdir"
  -B "${PROJECT_ROOT}:/base"
  -B "${LICENSES_DIR}:/opts"
  -B "${scratchdir}:/scratch"
  "$FMRIPREP_IMAGE"
  /base/bids /base/derivatives/fmriprep
  participant --participant_label "$sub"
  --stop-on-first-crash
  --me-output-echos
  --output-spaces "${output_spaces[@]}"
  --cifti-output "$FMRIPREP_CIFTI_DENSITY"
  --bids-filter-file /base/code/fmriprep_config.json
  --skip-bids-validation
  --fs-license-file /opts/fs_license.txt
  --fs-subjects-dir /base/derivatives/freesurfer
  -w /scratch
)

printf 'fMRIPrep command:'
printf ' %q' "${cmd[@]}"
printf '\n'
if ((dry_run)); then
  echo "Dry run: not launching fMRIPrep."
  python3 "${scriptdir}/check_pipeline_state.py" fmriprep-complete "$bidsdir" "$derivdir" "$sub" --list || true
  exit 0
fi

mkdir -p "$derivdir" "$freesurferdir" "$scratchdir"
rf1_require_file "$FMRIPREP_IMAGE"
rf1_require_file "${LICENSES_DIR}/fs_license.txt"
"${cmd[@]}"
python3 "${scriptdir}/check_pipeline_state.py" fmriprep-complete "$bidsdir" "$derivdir" "$sub"
