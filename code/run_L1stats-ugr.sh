#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"
nruns=2
task=ugr # edit if necessary

for ppi in 0 aIns; do #0
	 #for sub in `cat ${scriptdir}/sublist_DD128.txt`; do
	 for sub in 10636; do	
	#for sub in `cat ${scriptdir}/sublist_mkExist.txt`; do	
	  for run in `seq $nruns`; do

	  	# Manages the number of jobs and cores
	  	SCRIPTNAME=${scriptdir}/L1stats-ugr.sh
	  	NCORES=50
	  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
	    		sleep 5s
	  	done
	  	bash $SCRIPTNAME $sub $run $ppi &
			sleep 1s
	  done
	done
done
