#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
if [[ -z "${FSLDIR:-}" ]] || ! command -v melodic >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v melodic >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOCK="$PWD/code/insula_ica/output/unsmoothed-group-ica-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-unsmoothed-group-ica.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/run_unsmoothed_ica.py --jobs 64 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started corrected unsmoothed insula ICA: PID $!"
echo "64 preparation slots, 18 fixed fits; six automatic fits run serially with 8 threads and 32 GiB RAM reservations."
echo "Both residual scales; left/right/bilateral; dimensions 3/5/7/automatic. Native MELODIC variance normalization enabled."
echo "TR-preserving, physically masked 4D inputs. No fMRIPrep, FEAT or FILM reruns."
echo "Results/failures publish to GitHub. Log: $INSULA_LOG"
echo "Status: $INSULA_PYTHON code/insula_ica/status_unsmoothed_ica.py --brief"
