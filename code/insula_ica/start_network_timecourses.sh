#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
export FSLDIR="${FSLDIR:-/usr/local/fsl}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$INSULA_PYTHON" code/insula_ica/remove_withdrawn_ica_projections.py
INSULA_LOG=$(mktemp /tmp/rf1-network-timecourses.XXXXXX)
nohup "$INSULA_PYTHON" -u code/insula_ica/run_network_timecourse_analysis.py \
    --jobs 8 --publish "$@" > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started standardized network OLS task bars: PID $!"
echo "Three-run pilot (Trust, UGR, corrected TR), then all 490 runs. Log: $INSULA_LOG"
echo "Reuses all six primary ICA bases. No FILM, connectivity, randomise or preprocessing jobs."
echo "Status: $INSULA_PYTHON code/insula_ica/status_network_timecourses.py --brief"
