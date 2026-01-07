#!/bin/bash

# ensure paths are correct irrespective from where user runs the script
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
basedir="$(dirname "$scriptdir")"


nruns=2
task=trust # edit if necessary
ppi=0 # could also be dACC or AI, just check paths and pre-registration

for sub in `cat ${scriptdir}/sublist_DD128.txt`; do	
	  for run in `seq $nruns`; do

	  	# Manages the number of jobs and cores
	  	SCRIPTNAME=${scriptdir}/L1stats-trust.sh
	  	NCORES=50
	  	while [ $(ps -ef | grep -v grep | grep $SCRIPTNAME | wc -l) -ge $NCORES ]; do
	    		sleep 5s
	  	done
	  	bash $SCRIPTNAME $sub $run $ppi &
		sleep 1s
		
	  done
	done
done
