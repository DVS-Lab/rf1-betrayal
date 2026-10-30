#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
PPI_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$PPI_PYTHON"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p code/ppi_seed_audit/output
PPI_LOG=$(mktemp /tmp/rf1-PPI-seed-audit.XXXXXX)
printf '%s\n' "$PPI_LOG" > code/ppi_seed_audit/output/launcher.logpath
nohup "$PPI_PYTHON" -u code/ppi_seed_audit/audit.py --publish --max-bolds 4 "$@" > "$PPI_LOG" 2>&1 < /dev/null &
echo "Original PPI seed audit started: PID $! | log: $PPI_LOG"
echo "Reads original scripts/designs/physiological files; no PPI, FEAT or preprocessing rerun."
echo "If original referenced BOLD survives, compares masks in at most four serial 4D reads; no replacement BOLD."
echo "Results and small source archives publish automatically to insula-ica."
echo "Status: $PPI_PYTHON code/ppi_seed_audit/audit.py --brief"
