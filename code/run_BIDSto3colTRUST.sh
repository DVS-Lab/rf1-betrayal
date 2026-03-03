#!/usr/bin/env bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

baseout=${maindir}/derivatives/fsl/EVfiles
mkdir -p ${baseout}

# Loop over subjects in bids
for subdir in ${maindir}/bids/sub-*; do

    sub=$(basename "$subdir")
    sub=${sub#sub-}

    echo "Processing subject $sub"

    for run in 1 2; do

        input=${maindir}/bids/sub-${sub}/func/sub-${sub}_task-trust_run-${run}_events.tsv

        if [ -f "$input" ]; then
            outdir=${baseout}/sub-${sub}/trust
            mkdir -p "$outdir"

            outprefix=${outdir}/run-${run}

            bash ${scriptdir}/BIDSto3colTRUST.sh "$input" "$outprefix"
        else
            echo "Missing file: $input"
        fi

    done

done
