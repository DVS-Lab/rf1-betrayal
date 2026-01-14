#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# the "type" variable below is setting a path inside the main script
type=act
for sub in `cat ${scriptdir}/sublist_DD128.txt`; do
  	SCRIPTNAME=${maindir}/code/L2stats-ugr.sh
  	NCORES=50
  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
    		sleep 1s
  	done
  	bash $SCRIPTNAME $sub 01 $type &
  	sleep 1s
done
