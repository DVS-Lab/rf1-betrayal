#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
INSULA_SOURCE=${1:-code/insula_ica/records/linux2-trust-stock-filter-20261028T165853Z-314305}
INSULA_JOBS=${INSULA_JOBS:-64}
test -x "$INSULA_PYTHON"
test -f "$INSULA_SOURCE/state.json"
if [[ -z "${FSLDIR:-}" ]] || ! command -v feat >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v feat >/dev/null
INSULA_LOCK="$PWD/code/insula_ica/output/trust-stock-filter-resume-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Resume launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-trust-stock-resume.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_SOURCE INSULA_LOCK INSULA_JOBS
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/resume_trust_stock_filter.py \
        --record "$INSULA_SOURCE" --jobs "$INSULA_JOBS" --fmriprep-jobs 2 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started resume/handoff: PID $!"
echo "Progress/errors: tail -n 15 '$INSULA_LOG'"
echo "$INSULA_JOBS worker slots; existing FEAT jobs retained; completed fMRIPrep/FEAT/QC reused."
echo "New fMRIPrep launches capped at two, with available-memory checks for new stages."
