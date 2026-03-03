#!/usr/bin/env bash

# Ensure paths are relative to rf1-betrayal
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

datadir=${maindir}/bids
baseout=${maindir}/derivatives/fsl/EVfiles

mkdir -p ${baseout}

sub=$1

for run in 1 2; do

  input=${datadir}/sub-${sub}/func/sub-${sub}_task-ugr_run-${run}_events.tsv
  output=${baseout}/sub-${sub}/ugr
  mkdir -p ${output}

  if [ -f "$input" ]; then

    # Call local BIDSto3col.sh
    bash ${scriptdir}/BIDSto3colUGR.sh -h "Offer" "$input" "${output}/run-${run}"

    # Rename _pmod to _constant
    for file in ${output}/run-${run}_*_pmod.txt; do
      [ -e "$file" ] && mv "$file" "${file/_pmod/_constant}"
    done

  else
    echo "Missing file: $input"
  fi

done
