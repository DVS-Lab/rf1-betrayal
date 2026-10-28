#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
if [[ -z "${FSLDIR:-}" ]] || ! command -v feat >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v feat >/dev/null
INSULA_LOCK="$PWD/code/insula_ica/output/registration-sub-10834-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-registration-sub10834.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/repair_trust_registration.py --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started sub-10834 correction candidate: PID $!"
echo "Progress/errors: tail -n 12 '$INSULA_LOG'"
echo "One participant, both Trust runs; original confounds preserved. Review records push at completion."
