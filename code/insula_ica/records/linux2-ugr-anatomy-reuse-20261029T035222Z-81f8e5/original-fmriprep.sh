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
lockroot="${PROJECT_ROOT}/logs/locks"
lockdir="${lockroot}/fmriprep-sub-${sub}.lock"
scratchdir="${SCRATCH_ROOT}/$(whoami)/fmriprep-sub-${sub}"
if [[ "$sub" == "10668" && -f "${bidsdir}/sub-10668/ses-01/.rf1-10668-repair.json" ]]; then
  scratchdir="${scratchdir}-sharedreward-repair-v1"
fi
fmriprep_nprocs="${FMRIPREP_NPROCS:-$FMRIPREP_TOTAL_NPROCS}"
fmriprep_omp_nthreads="${FMRIPREP_OMP_NTHREADS:-8}"
fmriprep_mem_mb="${FMRIPREP_MEM_MB:-$FMRIPREP_TOTAL_MEM_MB}"
IFS=' ' read -r -a output_spaces <<< "$FMRIPREP_OUTPUT_SPACES"
lock_acquired=0

bids_subject_dir="${bidsdir}/sub-${sub}"
if [[ ! -d "$bids_subject_dir" ]]; then
  if ((dry_run)); then
    echo "SKIP BIDS subject not found: $bids_subject_dir"
    echo "Run prepdata before fMRIPrep."
    exit 0
  fi
  rf1_require_dir "$bids_subject_dir"
fi

if [[ "$overwrite" -ne 1 ]] && python3 "${scriptdir}/check_pipeline_state.py" fmriprep-complete "$bidsdir" "$derivdir" "$sub" >/dev/null; then
  echo "sub-${sub} already has practical fMRIPrep completion outputs; skipping"
  exit 0
fi

mkdir -p "$lockroot"
if ! mkdir "$lockdir" 2>/dev/null; then
  lock_pid="$(cat "${lockdir}/pid" 2>/dev/null || true)"
  lock_host="$(cat "${lockdir}/host" 2>/dev/null || true)"
  this_host="$(hostname 2>/dev/null || true)"
  if [[ -n "$lock_pid" && "$lock_host" == "$this_host" ]] && ! kill -0 "$lock_pid" 2>/dev/null; then
    echo "Removing stale fMRIPrep lock for sub-${sub}: $lockdir"
    rf1_remove_tree_under "$lockroot" "$lockdir"
    mkdir "$lockdir"
  else
    echo "LOCKED (skipping): fMRIPrep appears active for sub-${sub}: $lockdir"
    exit 0
  fi
fi
lock_acquired=1
{
  echo "$$" > "${lockdir}/pid"
  hostname > "${lockdir}/host" 2>/dev/null || true
  date -Is > "${lockdir}/started"
}
cleanup_lock() {
  if ((lock_acquired)); then
    rf1_remove_tree_under "$lockroot" "$lockdir" || true
  fi
}
trap cleanup_lock EXIT

export APPTAINERENV_TEMPLATEFLOW_HOME=/opt/templateflow
export APPTAINERENV_MPLCONFIGDIR=/opt/mplconfigdir

cmd=(
  apptainer run --cleanenv
  -B "${TEMPLATEFLOW_HOME}:/opt/templateflow"
  -B "${MPLCONFIGDIR_HOST}:/opt/mplconfigdir"
  -B "${PROJECT_ROOT}:/base"
  -B "${LICENSES_DIR}:/opts"
  -B "${scratchdir}:/scratch"
  "$FMRIPREP_IMAGE"
  /base/bids /base/derivatives/fmriprep
  participant --participant_label "$sub"
  --stop-on-first-crash
  --nprocs "$fmriprep_nprocs"
  --omp-nthreads "$fmriprep_omp_nthreads"
  --mem "$fmriprep_mem_mb"
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
subject_log="${PROJECT_ROOT}/logs/runs/$(date +%Y%m%d-%H%M%S)_fmriprep-sub-${sub}.log"
mkdir -p "$(dirname "$subject_log")"
echo "fMRIPrep subject log: $subject_log"
set +e
"${cmd[@]}" 2>&1 | tee "$subject_log"
fmriprep_status=${PIPESTATUS[0]}
set -e
if ((fmriprep_status != 0)); then
  echo "fMRIPrep container failed for sub-${sub} (exit ${fmriprep_status}); see $subject_log" >&2
  exit "$fmriprep_status"
fi
python3 "${scriptdir}/check_pipeline_state.py" fmriprep-complete "$bidsdir" "$derivdir" "$sub"
