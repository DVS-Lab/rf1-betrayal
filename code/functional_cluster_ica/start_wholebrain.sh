#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
test "$(git branch --show-current)" = insula-ica
NETWORK_PYTHON=/ZPOOL/data/tools/anaconda/tug87422/envs/tedana-26.0.3/bin/python
test -x "$NETWORK_PYTHON"
export FSLDIR="${FSLDIR:-/usr/local/fsl}"
if ! command -v randomise >/dev/null 2>&1; then
    set +u
    source "$FSLDIR/etc/fslconf/fsl.sh"
    set -u
    export PATH="$FSLDIR/bin:$PATH"
fi
command -v randomise >/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/rf1-functional-network-mpl}
mkdir -p code/insula_ica/output "$MPLCONFIGDIR"
NETWORK_LOG=$(mktemp /tmp/rf1-functional-wholebrain.XXXXXX)
printf '%s\n' "$NETWORK_LOG" > code/insula_ica/output/functional-wholebrain-launcher.logpath
nohup "$NETWORK_PYTHON" -u code/functional_cluster_ica/wholebrain.py \
    --jobs 32 --randomise-jobs 32 --permutations 5000 --publish "$@" \
    > "$NETWORK_LOG" 2>&1 < /dev/null &
echo "Whole-brain dual regression launched with nohup: PID $!"
echo "Log: $NETWORK_LOG"
echo "32 RAM-gated run workers, then 32 single-thread randomise jobs (5000 permutations)."
echo "3/5 components, sm0/sm5, both existing residual scales: eight bases, 32 components."
echo "Reuses completed ICA/FEAT inputs and all verified network checkpoints on resume."
echo "Anatomy, paired task contrasts, group maps and TFCE results publish to GitHub."
echo "Status: $NETWORK_PYTHON code/functional_cluster_ica/status_wholebrain.py --brief"
