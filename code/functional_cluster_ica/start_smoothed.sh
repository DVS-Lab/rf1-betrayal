#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
MATCHED_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$MATCHED_PYTHON"
export FSLDIR="${FSLDIR:-/usr/local/fsl}"
if ! command -v melodic >/dev/null 2>&1; then
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v feat >/dev/null
command -v melodic >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/rf1-functional-cluster-mpl}
mkdir -p code/insula_ica/output "$MPLCONFIGDIR"
MATCHED_LOG=$(mktemp /tmp/rf1-functional-sm5.XXXXXX)
printf '%s\n' "$MATCHED_LOG" > code/insula_ica/output/functional-sm5-launcher.logpath
nohup "$MATCHED_PYTHON" -u code/functional_cluster_ica/smoothed.py --jobs 64 --publish "$@" \
    > "$MATCHED_LOG" 2>&1 < /dev/null &
echo "Matched smoothed-arm submitter launched with nohup: PID $!"
echo "Log: $MATCHED_LOG"
echo "Waits for the running sm-0 arm; no interruption or unsmoothed refit."
echo "64 RAM-gated FEAT/preparation slots; six fixed ICA slots; automatic PCA serialized."
echo "Reuses verified matching sm-5 FEATs; only missing/mismatched fits get separate outputs."
echo "Both residual scales × dimensions 3/5/7/0; same mask, raw task cache, five-bar model and mean references."
echo "Results, paired comparisons and failures publish to GitHub."
echo "Status: $MATCHED_PYTHON code/functional_cluster_ica/status_smoothed.py --brief"
