#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
if [[ -z "${FSLDIR:-}" ]] || ! command -v feat >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v feat >/dev/null
command -v apptainer >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOCK="$PWD/code/insula_ica/output/ugr-anatomy-reuse-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect its pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-ugr-anatomy-reuse.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/recover_ugr_unsmoothed.py --fmriprep-jobs 32 --feat-jobs 64 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started automated UGR pilot → expansion → missing sm-0 FEAT recovery: PID $!"
echo "Queue limits: 32 fMRIPrep subjects (8 threads each), 64 FEAT workers; shared CPU/RAM gate."
echo "Pilot checks anatomy reuse, BOLD/reference/coverage and the unchanged joint FEAT model before expansion."
echo "Existing Trust sm-0 fits are retained. Records publish after the pilot and at completion/failure."
echo "Log: $INSULA_LOG"
echo "Compact progress: $INSULA_PYTHON code/insula_ica/status_ugr_recovery.py --brief"
