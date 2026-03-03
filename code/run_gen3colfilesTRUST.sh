#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"

for sub in `ls -d ${basedir}/bids/sub-*/`; do

          sub=${sub#*sub-}
          sub=${sub%/}
	bash ${scriptdir}/gen3colfiles.sh $sub
done
