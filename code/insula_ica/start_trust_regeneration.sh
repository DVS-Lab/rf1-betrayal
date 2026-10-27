#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
INSULA_AUDIT="$PWD/code/insula_ica/output/linux2-audit-20261026T213907Z-dc0e41/audit"
INSULA_OUTPUT="$PWD/derivatives/insula_ica/trust-rebuild-24-linux2-audit-20261026T213907Z-dc0e41"
test -x "$INSULA_PYTHON"
test -f "$INSULA_AUDIT/manifest.tsv"
test -f /ZPOOL/data/tools/fmriprep-24.1.1.simg
if [[ -z "${FSLDIR:-}" ]] || ! command -v feat >/dev/null 2>&1; then
    if [[ -r /usr/local/fsl/etc/fslconf/fsl.sh ]]; then
        export FSLDIR=/usr/local/fsl
        set +u
        source "$FSLDIR/etc/fslconf/fsl.sh"
        set -u
        export PATH="$FSLDIR/bin:$PATH"
    else
        echo "FSL setup is missing; source the normal linux2 FSL setup first."
        exit 1
    fi
fi
command -v feat >/dev/null
test -n "${FSLDIR:-}"
INSULA_LOCK="$PWD/code/insula_ica/output/trust-regeneration-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-trust-regeneration.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_AUDIT INSULA_OUTPUT INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/regenerate_trust.py \
        --audit-dir "$INSULA_AUDIT" --output-root "$INSULA_OUTPUT" --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started Trust regeneration: PID $!"
echo "Progress/errors: tail -n 20 '$INSULA_LOG'"
echo "Eight concurrent subjects; completed outputs reused. Records commit/push at completion."
