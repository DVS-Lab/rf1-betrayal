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
command -v feat_model >/dev/null
command -v apptainer >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOCK="$PWD/code/insula_ica/output/ugr-timing-repair-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-ugr-timing-repair.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/repair_ugr_timing.py
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started corrected-TR UGR rerun: PID $!"
echo "Audits all 490 run timings and both Trust runs for sub-10817."
echo "Reuses anatomy; fits corrected sm-0 and paired sm-5 TR controls; reuses 489 finished FEATs."
echo "Results and errors publish to insula-ica. Log: $INSULA_LOG"
echo "Progress: $INSULA_PYTHON code/insula_ica/status_ugr_timing_repair.py --brief"
