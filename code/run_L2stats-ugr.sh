#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
maindir="$(dirname "$scriptdir")"

# the "type" variable below is setting a path inside the main script
for type in "ppi_seed-pTPJ"; do # act ppi_seed-pTPJ
	#for sub in 10700; do
	for sub in `cat ${scriptdir}/sublist_sans.txt`; do
		# Manage the number of jobs and cores # `cat ${scriptdir}/sublist_all.txt
  	SCRIPTNAME=${maindir}/code/L2stats.sh
  	NCORES=20
  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
    		sleep 1s
  	done
  	bash $SCRIPTNAME $sub $type &
  	sleep 1s

	done
done
