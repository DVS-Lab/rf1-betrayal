#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
UGR_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$UGR_PYTHON"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
UGR_AUDIT_LOG=$(mktemp /tmp/rf1-original-group-audit.XXXXXX)
nohup "$UGR_PYTHON" -u code/ugr_reanalysis/audit_inputs.py --publish > "$UGR_AUDIT_LOG" 2>&1 < /dev/null &
echo "Read-only group-map / original UGR design audit started: PID $!"
echo "Log: $UGR_AUDIT_LOG"
echo "Results publish automatically to insula-ica. No FEAT or group rerun launched."
