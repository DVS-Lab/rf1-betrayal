#!/usr/bin/env bash
# Launch quietly, saving and publishing the diagnostic record on completion.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
INSULA_REBUILD="$PWD/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41"
test -x "$INSULA_PYTHON"
test -f "$INSULA_REBUILD/plan-pilot.json"
INSULA_LOCK="$PWD/code/insula_ica/output/bold-controls-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "A control launcher lock already exists: $INSULA_LOCK"
    echo "Inspect its pid/log before starting another experiment."
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-bold-controls.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_REBUILD INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    "$INSULA_PYTHON" -u code/insula_ica/run_bold_controls.py \
        --rebuild-root "$INSULA_REBUILD" --jobs 8 --execute --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
INSULA_PID=$!
printf '%s\n' "$INSULA_PID" > "$INSULA_LOCK/pid"
echo "Started controls: PID $INSULA_PID"
echo "Progress/errors: tail -n 20 '$INSULA_LOG'"
echo "Full results will be committed and pushed to insula-ica on completion."
