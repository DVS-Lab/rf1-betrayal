#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
FUNCTIONAL_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$FUNCTIONAL_PYTHON"
export FSLDIR="${FSLDIR:-/usr/local/fsl}"
if ! command -v melodic >/dev/null 2>&1; then
    test -f "$FSLDIR/etc/fslconf/fsl.sh"
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v melodic >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/rf1-functional-cluster-mpl}
mkdir -p code/insula_ica/output "$MPLCONFIGDIR"
FUNCTIONAL_LOG=$(mktemp /tmp/rf1-functional-cluster.XXXXXX)
printf '%s\n' "$FUNCTIONAL_LOG" > code/insula_ica/output/functional-cluster-launcher.logpath
nohup "$FUNCTIONAL_PYTHON" -u code/functional_cluster_ica/run.py --jobs 32 --ica-jobs 6 --publish "$@" \
  > "$FUNCTIONAL_LOG" 2>&1 < /dev/null &
echo "Functional-cluster workflow launched with nohup: PID $!"
echo "Log: $FUNCTIONAL_LOG"
echo "32 RAM-gated preparation/task workers; six fixed ICA slots; automatic PCA serialized."
echo "Mean-augmented insula plots publish first; functional ICA + five-bar task plots follow."
echo "Completed stages resume automatically; no anatomical ICA, fMRIPrep, FEAT or FILM reruns."
echo "Status: $FUNCTIONAL_PYTHON code/functional_cluster_ica/status.py --brief"
