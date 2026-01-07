#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"


for sub in `cat ${scriptdir}/sublist_DD128.txt`; do
	# Manage the number of jobs and cores # `cat ${scriptdir}/sublist_all.txt
  	SCRIPTNAME=${maindir}/code/L2stats-trust.sh
  	NCORES=80
  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
    		sleep 1s
  	done
  	bash $SCRIPTNAME $sub 01 $type &
  	sleep 1s

	done
done
