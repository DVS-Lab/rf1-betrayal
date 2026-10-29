#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
if [[ -z "${FSLDIR:-}" ]]; then export FSLDIR=/usr/local/fsl; fi
test -f "$FSLDIR/data/atlases/HarvardOxford-Cortical.xml"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOCK="$PWD/code/insula_ica/output/completed-ICA-exploration-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Exploration launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-ICA-exploration.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/explore_completed_ica.py --jobs 8 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started completed-ICA exploration: PID $!"
echo "Publishes completed fixed-fit maps first, then all primary task projections, then whole-brain background connectivity."
echo "Uses all six five-component bases; eight I/O workers with RAM gating. No ICA/FEAT/FILM jobs or existing inputs changed."
echo "Automatic ICA continues independently. Log: $INSULA_LOG"
echo "Status: $INSULA_PYTHON code/insula_ica/status_ica_exploration.py --brief"
