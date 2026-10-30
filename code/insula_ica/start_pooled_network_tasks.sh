#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOG=$(mktemp /tmp/rf1-pooled-network-tasks.XXXXXX)
nohup "$INSULA_PYTHON" -u code/insula_ica/pooled_network_tasks.py --jobs 16 --publish \
    > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started pooled UGR task model: PID $! | log: $INSULA_LOG"
echo "Uses cached timecourses; refits 249 small UGR models x six bases, reuses all Trust estimates."
echo "Publishes all five ICs x five bars, all ROIs/bases. No 4D reads, exclusions, FEAT, ICA or randomise."
