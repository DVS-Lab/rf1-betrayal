#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
UGR_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$UGR_PYTHON"
export FSLDIR=/usr/local/fsl
set +u
source "$FSLDIR/etc/fslconf/fsl.sh"
set -u
export PATH="$FSLDIR/bin:$PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p code/ugr_reanalysis/output
UGR_LOG=$(mktemp /tmp/rf1-ugr-activation-repair.XXXXXX)
printf '%s\n' "$UGR_LOG" > code/ugr_reanalysis/output/activation-launcher.logpath
nohup "$UGR_PYTHON" -u code/ugr_reanalysis/rerun_activation.py --jobs 8 --publish "$@" > "$UGR_LOG" 2>&1 < /dev/null &
echo "UGR activation repair started: PID $! | log: $UGR_LOG"
echo "One original-model group pilot, then remaining five groups concurrently; completed groups reused."
echo "Outputs: derivatives/ugr_reanalysis/<timestamp>/activation/<model>/"
echo "Status: $UGR_PYTHON code/ugr_reanalysis/rerun_activation.py --brief"
