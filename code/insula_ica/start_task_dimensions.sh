#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOG=$(mktemp /tmp/rf1-task-dimensions.XXXXXX)
mkdir -p code/insula_ica/output
printf '%s\n' "$INSULA_LOG" > code/insula_ica/output/task-dimensions-launcher.logpath
nohup "$INSULA_PYTHON" -u code/insula_ica/task_dimension_sensitivity.py --jobs 8 --publish "$@" \
    > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started 3/7-component task sensitivity plus insula means: PID $!"
echo "Log: $INSULA_LOG"
echo "Reuses finished ICA and 5-component task results; auto-resumes matching extraction/model checkpoints."
echo "One 4D read per unfinished run; eight RAM-gated readers. No ICA, FEAT, connectivity or randomise."
echo "Status: $INSULA_PYTHON code/insula_ica/status_task_dimensions.py --brief"
