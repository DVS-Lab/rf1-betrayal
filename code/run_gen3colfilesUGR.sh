#!/bin/bash

scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

for sub in `cat ${scriptdir}/sublist_n139.txt`; do
    bash ${scriptdir}/gen3colfilesUGR.sh $sub
done
