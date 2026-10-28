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
command -v film_gls >/dev/null
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
    "$INSULA_PYTHON" -u code/insula_ica/run_group_ica.py --jobs 64 --ica-jobs 18 --auto-ica-jobs 6 --publish
' > "$INSULA_LOG" 2>&1 < /dev/null &
echo "Started paired insula ICA pipeline: PID $!"
echo "Progress/errors: tail -n 15 '$INSULA_LOG'"
echo "N132; Trust241 + UGR249; all three reviewed runs retained."
echo "Five scale-validation pilots, then 64 preparation slots with RAM checks; reuse/recover native FILM residuals as needed."
echo "Left/right/bilateral × 3/5/7/auto(0) components × joint-BOLD/native-FILM: 24 fits."
echo "18 fixed-dimensionality fits up to 18 concurrent; then six auto fits with 16 GiB RAM reservations plus 8 GiB headroom (four threads per fit)."
