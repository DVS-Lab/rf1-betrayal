#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
INSULA_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$INSULA_PYTHON"
if [[ -z "${FSLDIR:-}" ]] || ! command -v melodic >/dev/null 2>&1; then
    export FSLDIR=/usr/local/fsl
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v melodic >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
INSULA_LOCK="$PWD/code/insula_ica/output/paired-group-ica-launch.lock"
mkdir -p "$(dirname "$INSULA_LOCK")"
if ! mkdir "$INSULA_LOCK" 2>/dev/null; then
    echo "Launch lock exists; inspect pid/log: $INSULA_LOCK"
    exit 1
fi
INSULA_LOG=$(mktemp /tmp/rf1-paired-group-ica.XXXXXX)
printf '%s\n' "$INSULA_LOG" > "$INSULA_LOCK/log"
export INSULA_PYTHON INSULA_LOCK
nohup bash -c '
    trap '\''rm -f "$INSULA_LOCK/pid" "$INSULA_LOCK/log"; rmdir "$INSULA_LOCK"'\'' EXIT
    printf "%s\n" "$$" > "$INSULA_LOCK/pid"
    "$INSULA_PYTHON" -u code/insula_ica/resume_group_ica.py --jobs 11 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started ICA-only resume: PID $!"
echo "Progress/errors: tail -n 15 '$INSULA_LOG'"
echo "Reuses verified 490-run preparation, 1584 masked 4D inputs, and all completed ICA fits."
echo "Only unfinished fits run: explicit seeds 20261028–20261032, up to 10000 iterations, unchanged concat/pow3/epsilon."
echo "11 queue slots; fixed fits reserve 2 GiB, auto fits 16 GiB, plus 8 GiB headroom. No FILM or input rewriting."
echo "Then a separate native-MELODIC variance-normalization pilot: right insula, five components, both residual types; uniformly rescaled copies only."
